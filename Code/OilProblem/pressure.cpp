//
// Created by Stefano Weidmann on 27.03.18.
//

#include "oilproblem.hpp"
#include "pressureInternal.hpp"

#ifdef __GNUC__
    #define likely(x)       __builtin_expect((x),1)
    #define unlikely(x)     __builtin_expect((x),0)
#else
    #define likely(x) x
    #define unlikely(x) x
#endif



Real computeTransmissibility(ConstMatrixRef lambdas, const CellIndex& fromCell, const CellIndex& toCell) {
    // transmissibility is harmonic mean of lambda coefficients

    const double lambdaFrom = fromCell(lambdas);
    const double lambdaTo = toCell(lambdas);

    return 2.0 * lambdaFrom * lambdaTo / (lambdaFrom + lambdaTo);
}



CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const int numberOfRows) {

    return {fromCell.linearIndex(numberOfRows), toCell.linearIndex(numberOfRows)};
}



SparseMatrix assembleTransmissibilityMatrix(ConstMatrixRef lambdas) {

    const int numberOfRows = int(lambdas.rows());
    const int numberOfCols = int(lambdas.cols());


    const int numberOfPairs = numberOfCols * numberOfCols;
    SparseMatrix transmissibilities(numberOfPairs, numberOfPairs);



    transmissibilities.reserve(Eigen::VectorXi::Constant(numberOfPairs, 5));

    CellIndex myself = {0, 0};

    const CellIndex wellCell = {0, numberOfCols-1};



    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {


            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows);
            constexpr static std::array<CellIndex::Direction, 2> directionsToCheck = {
                  CellIndex::Direction::SOUTH, CellIndex::Direction::WEST
            };

            for (const CellIndex::Direction direction: directionsToCheck) {
                if (unlikely(!myself.hasNeighbor(direction, numberOfRows, numberOfCols))) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);

                const CellIndex neighborToThemselves = pressureToTransmissibilityIndex(neighbor, neighbor, numberOfRows);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);
                const CellIndex neighborToMe = meToNeighbor.transpose();

                const Real currentTransmissibility = computeTransmissibility(lambdas, myself, neighbor);

                const bool iAmTheCellInTheWell = myself == wellCell;

                if (unlikely(iAmTheCellInTheWell)) {
                    // pressure of cell in the well is known
                    // By the book, it would be nicer to exclude the row and column of the cell in the well
                    // and add some terms to the rhs of the linear system of equations in terms of execution speed.
                    // This however is more implementation effort and requires returning the correction of the rhs,
                    // which is not nice.
                    // May change that later.
                    meToMyself(transmissibilities) = 1;
                } else {
                    meToMyself(transmissibilities) += currentTransmissibility;
                    meToNeighbor(transmissibilities) -= currentTransmissibility;
                }

                const bool neighborIsCellInTheWell = neighbor == wellCell;

                if (likely(!neighborIsCellInTheWell)) {
                    neighborToMe(transmissibilities) -= currentTransmissibility;
                    neighborToThemselves(transmissibilities) += currentTransmissibility;
                }
            }
        }
    }

    transmissibilities.makeCompressed();
    return transmissibilities;
}

VectorToBeMappedAsMatrix solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstMatrixRef sources) {
    const Eigen::Map<const Vector> sourcesAsVector(sources.data(), sources.size());

    SparseSolver solver;
    solver.compute(transmissibilities);

    Vector pressuresAsVector = solver.solve(sourcesAsVector);

    // Eigen map doesn't take ownership of the vector it mapped, so we have to pass it unfortunately
    return VectorToBeMappedAsMatrix(std::move(pressuresAsVector), int(sources.rows()), int(sources.cols()));
}

