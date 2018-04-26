//
// Created by Stefano Weidmann on 17.04.18.
//

#include "oilProblem.hpp"

[[deprecated]]
SparseMatrix makeMatrixColStochastic(SparseMatrix matrix) {
    for (int col = 0; col < matrix.cols(); ++col) {
        matrix.col(col) /= matrix.col(col).sum();
    }

    return matrix;
}


[[deprecated]]
Real getTransitionProbability(const int fromState, const int toState, const SparseMatrix& transitionProbabilitiesColStochastic) {
    return transitionProbabilitiesColStochastic.coeff(toState, fromState);
}


// TODO transition
void transitionState(RandomWalkState& currentState,
                                const SparseMatrix& pressureResidualsByPressures,
                                const SparseMatrix& pressureResidualsBySaturationsWater,
                                const SparseMatrix& saturationsWaterResidualsByPressure,
                                const SparseMatrix& saturationsWaterResidualsBySaturationsWater) {


    std::vector<Real> candidates(10);


    if (currentState.isAPressure) {
        // stay a pressure



    } else {

        ++currentState.currentTimelevel;
    }
}

AdjointState initialAdjointState(const int numberOfRows, const int numberOfCols, const int numberOfParameters, const BVectorSurrogate& b, const CMatrixSurrogate& c) {
    AdjointState initialState;

    const int numberOfRandomWalksPerPressureCell = 5;

    const int numberOfRandomWalksPerParameter = numberOfRandomWalksPerPressureCell;

    // We start randomWalksPerPressureCell random walks for every pressure state, which means
    // randomWalksPerPressureCell per cell
    // We hope that the saturation states are reached through the coupling of the equations

    initialState.estimators.resize(numberOfRandomWalksPerPressureCell, numberOfParameters);
    initialState.weights.resize(numberOfRandomWalksPerPressureCell, numberOfParameters);


    for (int parameterIndex = 0; parameterIndex <  numberOfParameters; ++parameterIndex) {
        const CellIndex cell = CellIndex::fromLinearIndex(parameterIndex, numberOfRows);

        std::vector<CellIndex> neighbors;

        constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
              CellIndex::Direction::EAST, CellIndex::Direction::WEST, CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH
        };

        for (const auto direction: directionsToCheck) {
            if (cell.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                neighbors.push_back(cell.neighbor(direction));
            }
        }

        for (int randomWalkIndex = 0; randomWalkIndex < numberOfRandomWalksPerPressureCell; ++randomWalkIndex) {
            const int whichCellIndex = randomWalkIndex % int(neighbors.size());
            const CellIndex neighborOrMyself = neighbors[whichCellIndex];

            const bool wantAPressure = true;

            initialState.weights(randomWalkIndex, parameterIndex) = numberOfRandomWalksPerPressureCell * c(neighborOrMyself, cell, wantAPressure);
            initialState.estimators(randomWalkIndex, parameterIndex) = initialState.weights(randomWalkIndex, parameterIndex) * b(neighborOrMyself, wantAPressure);
        }
    }

    return initialState;



}
