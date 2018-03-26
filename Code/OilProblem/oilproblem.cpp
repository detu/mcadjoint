//
// Created by Stefano Weidmann on 24.03.18.
//

#include "oilproblem.hpp"

#include <Eigen/IterativeLinearSolvers>

Real computeTransmissibility(ConstMatrixRef lambdas, const CellIndex& fromCell, const CellIndex& toCell) {
    // transmissibility is harmonic mean of lambda coefficients

    const double lambdaFrom = fromCell(lambdas);
    const double lambdaTo = toCell(lambdas);

    return 2.0 * lambdaFrom * lambdaTo / (lambdaFrom + lambdaTo);
}

CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const long numberOfRows) {

    return {fromCell.linearIndex(numberOfRows), toCell.linearIndex(numberOfRows)};
}



AssemblyOutput assembleTransmissibilityMatrix(ConstMatrixRef lambdas, const double pressureAtWell) {

    // TRANS MATRIX SINGULAR
    const int numberOfRows = int(lambdas.rows());
    const int numberOfCols = int(lambdas.cols());


    const long numberOfPairs = numberOfCols * numberOfCols;
    AssemblyOutput output {SparseMatrix(numberOfPairs-1, numberOfPairs-1), SparseVector(numberOfPairs-1)};
    SparseMatrix transmissibilities;



    transmissibilities.reserve(Eigen::VectorXi::Constant(numberOfPairs, 5));

    CellIndex currentPressureCell = {0, 0};

    const CellIndex wellCellIndex = {0, numberOfCols-1};



    for (currentPressureCell.j = 0; currentPressureCell.j < numberOfCols; ++currentPressureCell.j) {
        for (currentPressureCell.i = 0; currentPressureCell.i < numberOfRows; ++currentPressureCell.i) {


            const CellIndex meToMyself = pressureToTransmissibilityIndex(currentPressureCell, currentPressureCell, numberOfRows);
            constexpr static std::array<CellIndex::Direction, 2> directionsToCheck = {
                  CellIndex::Direction::SOUTH, CellIndex::Direction::WEST
            };

            for (const CellIndex::Direction direction: directionsToCheck) {
                if (!currentPressureCell.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const CellIndex neighbor = currentPressureCell.neighbor(direction);

                const CellIndex neighborToThemselves = pressureToTransmissibilityIndex(neighbor, neighbor, numberOfRows);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(currentPressureCell, neighbor, numberOfRows);
                const CellIndex neighborToMe = meToNeighbor.transpose();

                const Real currentTransmissibility = computeTransmissibility(lambdas, currentPressureCell, neighbor);
                neighborToThemselves(transmissibilities) += currentTransmissibility;
                meToNeighbor(transmissibilities) -= currentTransmissibility;

                const bool iAmTheCellAtTheWell = currentPressureCell == wellCellIndex;

                if (iAmTheCellAtTheWell) {
                    meToMyself(transmissibilities) = 1;
                }
                meToMyself(transmissibilities) += currentTransmissibility;
                neighborToMe(transmissibilities) -= currentTransmissibility;
            }
        }
    }

    transmissibilities.makeCompressed();
    return transmissibilities;
}

Matrix solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstMatrixRef sources) {
    const Eigen::Map<const Vector> sourcesAsVector(sources.data(), sources.size());

    Eigen::SimplicialLDLT<SparseMatrix> solver;
    solver.compute(transmissibilities);

    Vector pressuresAsVector = solver.solve(sourcesAsVector);

    return std::move(Eigen::Map<Matrix>(pressuresAsVector.data(), sources.rows(), sources.cols()));
}