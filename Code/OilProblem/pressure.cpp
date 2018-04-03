//
// Created by Stefano Weidmann on 27.03.18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>

#ifdef __GNUC__
    #define likely(x)       __builtin_expect((x),1)
    #define unlikely(x)     __builtin_expect((x),0)
#else
    #define likely(x) x
    #define unlikely(x) x
#endif



Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell) {
    // transmissibility is harmonic mean of total mobility coefficients

    const double lambdaFrom = fromCell(totalMobilities);
    const double lambdaTo = toCell(totalMobilities);

    return 2.0 * lambdaFrom * lambdaTo / (lambdaFrom + lambdaTo);
}

CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const int numberOfRows,
      const int numberOfCols) {


    const CellIndex wellCell = {0, numberOfCols-1};
    const int wellCellLinearIndex = wellCell.linearIndex(numberOfRows);

    return {fromCell.linearIndex(numberOfRows),
            toCell.linearIndex(numberOfRows)};
}



SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities) {

    const int numberOfRows = totalMobilities.rows();
    const int numberOfCols = totalMobilities.cols();


    const int numberOfPairs = numberOfCols * numberOfRows;
    SparseMatrix transmissibilities(numberOfPairs, numberOfPairs);


    transmissibilities.reserve(Eigen::VectorXi::Constant(transmissibilities.cols(), 6));

    CellIndex myself = {0, 0};

    const CellIndex wellCell = {0, numberOfCols-1};

    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {


            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows, numberOfCols);
            constexpr static std::array<CellIndex::Direction, 2> directionsToCheck = {
                  CellIndex::Direction::SOUTH, CellIndex::Direction::WEST
            };

            for (const CellIndex::Direction direction: directionsToCheck) {
                if (unlikely(!myself.hasNeighbor(direction, numberOfRows, numberOfCols))) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);

                const CellIndex neighborToThemselves = pressureToTransmissibilityIndex(neighbor, neighbor, numberOfRows,
                                                                                       numberOfCols);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows,
                                                                               numberOfCols);
                const CellIndex neighborToMe = meToNeighbor.transpose();

                const Real currentTransmissibility = computeTransmissibility(totalMobilities, myself, neighbor);

                const bool neighborIsCellInTheWell = wellCell == neighbor;
                const bool iAmCellInTheWell = wellCell == myself;

                // TODO SparseQr
                if (likely(!iAmCellInTheWell)) {
                    meToMyself(transmissibilities) += currentTransmissibility;

                    if (likely(!neighborIsCellInTheWell)) {
                        neighborToMe(transmissibilities) -= currentTransmissibility;
                        neighborToThemselves(transmissibilities) += currentTransmissibility;
                    }
                } else {
                    meToMyself(transmissibilities) = 1;
                }
            }
        }
    }

    transmissibilities.makeCompressed();
    return transmissibilities;
}


void solvePressurePoissonProblemInplace(const SparseMatrix& transmissibilities, ConstMatrixRef sources, const Real pressureAtWell, VectorRef result) {

    SparseSolver solver;
    solver.compute(transmissibilities);
    result = solver.solve(Eigen::Map<const Vector>(sources.data(), sources.size()));
    //std::cout << "Estimated error = " << solver.error() << "\n";

    result.array() += pressureAtWell;
}

Vector assemblePressureSourceVector(const Real pressureWellNow, const int numberOfRows, const int numberOfCols) {
    const CellIndex wellCell = {0, numberOfCols - 1};

    Vector pressureSourceVector(numberOfRows * numberOfCols);


}