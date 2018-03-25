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

SparseMatrix assembleTransmissibilityMatrix(ConstMatrixRef lambdas) {

    // TRANS MATRIX SINGULAR
    const long numberOfRows = lambdas.rows();
    const long numberOfCols = lambdas.cols();


    const long numberOfPairs = numberOfCols * numberOfCols;
    SparseMatrix transmissibilities(numberOfPairs, numberOfPairs);



    transmissibilities.reserve(Eigen::VectorXi::Constant(numberOfPairs, 5));

    CellIndex currentPressureCell = {0, 0};


    // TODO Make pressure matrix regular
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
                meToMyself(transmissibilities) += currentTransmissibility;
                neighborToThemselves(transmissibilities) += currentTransmissibility;
                meToNeighbor(transmissibilities) -= currentTransmissibility;
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