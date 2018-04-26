//
// Created by Stefano Weidmann on 17.04.18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>
#include <random>
#include <stefCommonHeaders/xoroshiro.h>

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
                                const SparseMatrix& saturationsWaterResidualsBySaturationsWater,
                                const int numberOfRows, const int numberOfCols, Rng& rng) {


    struct Candidate {
        const CellIndex cellIndex;
        const bool isAPressure;
    };

    std::vector<Real> candidateUnnormalizedProbabilities(10);
    std::vector<Candidate> candidates(10);

    constexpr bool iAmAPressure = true;
    constexpr bool iAmASaturation = !iAmAPressure;

    const SparseMatrix& pressureResidualsDerived = (currentState.isAPressure? pressureResidualsByPressures: pressureResidualsBySaturationsWater);
    const SparseMatrix& saturationWaterResidualsDerived = currentState.isAPressure? saturationsWaterResidualsByPressure: saturationsWaterResidualsBySaturationsWater;

    if (currentState.isAPressure) {
        // for pressures there's the probability to stay at the same place
        const CellIndex meToMyself = pressureToTransmissibilityIndex(currentState.cell, currentState.cell, numberOfRows);
        const Real unnormalizedProbabilityOfStayingHere = 1 - std::abs(meToMyself(pressureResidualsByPressures));
        ASSERT(unnormalizedProbabilityOfStayingHere > 0);
        candidates.push_back({currentState.cell, iAmAPressure});
        candidateUnnormalizedProbabilities.push_back(unnormalizedProbabilityOfStayingHere);
    }

    // add neighbors
    for (const auto neighbor: currentState.cell.neighbors(numberOfRows, numberOfCols)) {
        const CellIndex neighborToMe = pressureToTransmissibilityIndex(neighbor, currentState.cell, numberOfRows);
        const Real candidateUnnormalizedProbabilityPressure = std::abs(neighborToMe(pressureResidualsDerived));
        ASSERT(candidateUnnormalizedProbabilityPressure > 0);
        const Real candidateUnnormalizedProbabilitySaturationWater = std::abs(neighborToMe(saturationWaterResidualsDerived));
        ASSERT(candidateUnnormalizedProbabilitySaturationWater > 0);

        candidates.push_back({neighbor, iAmAPressure});
        candidates.push_back({neighbor, iAmASaturation});
        candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilityPressure);
        candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilitySaturationWater);

    }

    std::discrete_distribution<int> choose(candidateUnnormalizedProbabilities.cbegin(), candidateUnnormalizedProbabilities.cend());
    const int nextIndex = choose(rng);
    const Candidate& chosenCandidate = candidates[nextIndex];
    if (!(currentState.isAPressure && chosenCandidate.isAPressure)) {
        ++currentState.currentTimelevel;
    }

    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;
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
