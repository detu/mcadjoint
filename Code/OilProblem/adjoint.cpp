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
// Senke bei pressures drill cell!!
void transitionState(RandomWalkState& currentState, AdjointState& adjointState,
                                const SparseMatrix& pressureResidualsByPressures,
                                const SparseMatrix& pressureResidualsBySaturationsWater,
                                const SparseMatrix& saturationsWaterResidualsByPressure,
                                const SparseMatrix& saturationsWaterResidualsBySaturationsWater,
                                const int numberOfRows, const int numberOfCols, Rng& rng) {


    struct Candidate {
        const CellIndex cellIndex;
        const bool isAPressure;
        const Real correspondingEntryOfAMatrix;
        const Real correspondingEntryOfBVector;
    };

    std::vector<Real> candidateUnnormalizedProbabilities(10);
    std::vector<Candidate> candidates(10);

    constexpr bool iAmAPressure = true;
    constexpr bool iAmASaturation = !iAmAPressure;

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    const Real importanceOfDrillCellForPressure = 1000;

    const SparseMatrix& pressureResidualsDerived = (currentState.isAPressure? pressureResidualsByPressures: pressureResidualsBySaturationsWater);
    const SparseMatrix& saturationWaterResidualsDerived = currentState.isAPressure? saturationsWaterResidualsByPressure: saturationsWaterResidualsBySaturationsWater;

    int numberOfPlacesToGo = 0;
    if (currentState.isAPressure) {
        // for pressures there's the probability to stay at the same place
        const CellIndex meToMyself = pressureToTransmissibilityIndex(currentState.cell, currentState.cell, numberOfRows);
        const Real unnormalizedProbabilityOfStayingHere = 1;/*1 - std::abs(meToMyself(pressureResidualsByPressures))*/;
        ASSERT(unnormalizedProbabilityOfStayingHere > 0);
        candidates.push_back({currentState.cell, iAmAPressure, 1 - std::abs(meToMyself(pressureResidualsByPressures))});
        candidateUnnormalizedProbabilities.push_back(unnormalizedProbabilityOfStayingHere);
        ++numberOfPlacesToGo;
    }

    // add neighbors
    for (const auto neighbor: currentState.cell.neighbors(numberOfRows, numberOfCols)) {
        const CellIndex neighborToMe = pressureToTransmissibilityIndex(neighbor, currentState.cell, numberOfRows);
        const Real candidateUnnormalizedProbabilityPressure = neighbor == drillCell? importanceOfDrillCellForPressure: 1; /*std::abs(neighborToMe(pressureResidualsDerived))*/;
        ASSERT(candidateUnnormalizedProbabilityPressure > 0);
        const Real candidateUnnormalizedProbabilitySaturationWater = 1 /*std::abs(neighborToMe(saturationWaterResidualsDerived))*/;
        ASSERT(candidateUnnormalizedProbabilitySaturationWater > 0);

        candidates.push_back({neighbor, iAmAPressure, -neighborToMe(pressureResidualsDerived)});
        candidates.push_back({neighbor, iAmASaturation, -neighborToMe(saturationWaterResidualsDerived)});
        candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilityPressure);
        candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilitySaturationWater);
        numberOfPlacesToGo += 2;

    }

    std::discrete_distribution<int> choose(candidateUnnormalizedProbabilities.cbegin(), candidateUnnormalizedProbabilities.cend());
    const int nextIndex = choose(rng);
    const Candidate& chosenCandidate = candidates[nextIndex];
    if (!(currentState.isAPressure && chosenCandidate.isAPressure)) {
        ++currentState.currentTimelevel;
    }
    // update W (pg. 6199, top)
    if (chosenCandidate.isAPressure && chosenCandidate.)
    currentState.W *= numberOfPlacesToGo * chosenCandidate.transitionEntry;
    // update D (pg. 6199, top)
    currentState.D += currentState.W * chosenCandidate.correspondingEntryOfBVector;


    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;
}

std::vector<RandomWalkState> initializeRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters, const BVectorSurrogate& b, const CMatrixSurrogate& c) {

    const int numberOfRandomWalksPerPressureCell = 5;

    const int numberOfRandomWalksPerParameter = numberOfRandomWalksPerPressureCell;
    const int numberOfRandomWalks = numberOfRandomWalksPerParameter * numberOfParameters;

    // We start randomWalksPerPressureCell random walks for every pressure state, which means
    // randomWalksPerPressureCell per cell
    // We hope that the saturation states are reached through the coupling of the equations
    
    std::vector<RandomWalkState> randomWalks(numberOfRandomWalks);
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

        for (int localRandomWalkIndex = 0; localRandomWalkIndex < numberOfRandomWalksPerPressureCell; ++localRandomWalkIndex) {
            const int whichCellIndex = localRandomWalkIndex % int(neighbors.size());
            const CellIndex neighborOrMyself = neighbors[whichCellIndex];

            constexpr bool wantAPressure = true;
            RandomWalkState initialState;
            initialState.cell = neighborOrMyself;
            initialState.isAPressure = wantAPressure;
            initialState.currentTimelevel = 0;
            initialState.W = numberOfRandomWalksPerPressureCell * c(neighborOrMyself, cell, wantAPressure);
            initialState.D = initialState.weights(localRandomWalkIndex, parameterIndex) * b(neighborOrMyself, wantAPressure);

            randomWalks.push_back(std::move(initialState));
        }
    }

    return randomWalks;



}
