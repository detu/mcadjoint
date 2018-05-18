//
// Created by Stefano Weidmann on 17.05.18.
//

#include "regularization.hpp"
#include "darcyVelocity.hpp"

static constexpr Real regularizationWeight = 1e-3;

Real computeRegularizationPenalty(ConstMatrixRef logPermeabilities, const Real meshWidth) {
    return (computeXDerivative(logPermeabilities, meshWidth).cwiseAbs2().sum() + computeYDerivative(logPermeabilities, meshWidth).cwiseAbs2().sum()) * regularizationWeight;
}

Vector deriveRegularizationPenaltyByLogPermeabilities(ConstMatrixRef logPermeabilities, const Real meshWidth) {
    Vector derivatives = Vector::Zero(logPermeabilities.size());
    CellIndex cell = {0, 0};
    const int numberOfRows = logPermeabilities.rows();
    const int numberOfCols = logPermeabilities.cols();

    for (cell.j = 0; cell.j < numberOfCols; ++cell.j) {
        for (cell.i = 0; cell.i < numberOfRows; ++cell.i) {
            const int myLinearIndex = cell.linearIndex(numberOfRows);

            for (const CellIndex& neighbor: cell.neighbors(numberOfRows, numberOfCols)) {
                derivatives(myLinearIndex) += 2.0 * (cell(logPermeabilities) - neighbor(logPermeabilities)) / std::pow(meshWidth, 2) * regularizationWeight;
            }
        }
    }

    return derivatives;
}