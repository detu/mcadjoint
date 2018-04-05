//
// Created by Stefano Weidmann on 27.03.18.
//

#include "oilProblem.hpp"
#include "logging.hpp"
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
    LOGGER->debug("Starting to assemble system");

    const int numberOfRows = totalMobilities.rows();
    const int numberOfCols = totalMobilities.cols();


    const int numberOfPairs = numberOfCols * numberOfRows;
    SparseMatrix transmissibilities(numberOfPairs  + 1, numberOfPairs); // + 1 because we want to fix the pressure at the well to be zero


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


                if (wellCell != myself) {
                    meToMyself(transmissibilities) += currentTransmissibility;
                    neighborToMe(transmissibilities) -= currentTransmissibility;
                }

                if (wellCell != neighbor) {
                    meToNeighbor(transmissibilities) -= currentTransmissibility;
                    neighborToThemselves(transmissibilities) += currentTransmissibility;
                }
            }
        }
    }

    transmissibilities.coeffRef(transmissibilities.rows() - 1, wellCell.linearIndex(numberOfRows)) = 1;
    transmissibilities.makeCompressed();

    LOGGER->debug("Assembled system");
    return transmissibilities;
}


Vector projectSourcesIntoRange(ConstMatrixRef sources) {
    LOGGER->debug("Starting to project sources into range");
    Vector sourcesProjectedIntoRange(sources.size() + 1);
    sourcesProjectedIntoRange << Eigen::Map<const Vector>(sources.data(), sources.size()), 0;

    sourcesProjectedIntoRange.head(sources.size()).array() -= sourcesProjectedIntoRange.head(sources.size()).array().mean();

    LOGGER->debug("Projected sources into range");

    return sourcesProjectedIntoRange;
}

Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstVectorRef sourcesProjectedIntoRange, ConstVectorRef pressureGuess) {
    LOGGER->debug("Solving system");

    Eigen::LeastSquaresConjugateGradient<SparseMatrix> solver;
    //solver.setTolerance(1e-3);
    //solver.setMaxIterations(1);
    solver.compute(transmissibilities);

    const Vector result = solver.solveWithGuess(sourcesProjectedIntoRange, pressureGuess);

    LOGGER->debug("System solved");

    ASSERT(solver.info() == Eigen::Success);
    //std::cout << "Estimated error = " << solver.error() << "\n";

    return result;
}

Vector augmentSources(ConstMatrixRef sources) {

    Vector augmentedSources(sources.size() + 1);
    augmentedSources << Eigen::Map<const Vector>(sources.data(), sources.size()), 0;
    return augmentedSources;
}

#ifdef TODO
Vector assemblePressureSourceVector(const Real pressureWellNow, const int numberOfRows, const int numberOfCols) {
    const CellIndex wellCell = {0, numberOfCols - 1};

    Vector pressureSourceVector(numberOfRows * numberOfCols);


}
    #endif