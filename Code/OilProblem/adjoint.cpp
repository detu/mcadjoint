//
// Created by Stefano Weidmann on 17.04.18.
//

#include "adjoint.hpp"
#include "specialCells.hpp"
#include <stefCommonHeaders/assert.h>
#include <random>
#include <stefCommonHeaders/xoroshiro.h>
#include <iostream>
#include "logging.hpp"



bool transitionState(RandomWalkState& currentState, const BVectorSurrogate& b,
                                const SparseMatrix& pressureResidualsByPressures,
                                const SparseMatrix& pressureResidualsBySaturationsWater,
                                const SparseMatrix& saturationsWaterResidualsByPressure,
                                const SparseMatrix& saturationsWaterResidualsBySaturationsWater,
                                const int numberOfRows, const int numberOfCols, Rng& rng) {



    bool stillInTheSameTimestep;
    if (currentState.cell == CellIndex::invalidCell()) {
        stillInTheSameTimestep = false;
        return stillInTheSameTimestep;
    }

    ASSERT(currentState.cell.i < numberOfRows);
    ASSERT(currentState.cell.j < numberOfCols);
    ASSERT(currentState.cell.i >= 0);
    ASSERT(currentState.cell.j >= 0);


    struct Candidate {
        const CellIndex cellIndex;
        const bool isAPressure;
        const Real correspondingEntryOfAMatrix;
        const Real correspondingEntryOfBVector;
    };



    std::vector<Real> candidateUnnormalizedProbabilities;
    candidateUnnormalizedProbabilities.reserve(10);
    std::vector<Candidate> candidates;
    candidates.reserve(10);

    constexpr bool iAmAPressure = true;
    constexpr bool iAmASaturation = !iAmAPressure;
    constexpr bool neighborIsAPressure = iAmAPressure;
    constexpr bool neighborIsASaturation = iAmASaturation;

    const SparseMatrix& pressureResidualsDerived = (currentState.isAPressure? pressureResidualsByPressures: pressureResidualsBySaturationsWater);
    const SparseMatrix& saturationWaterResidualsDerived = currentState.isAPressure? saturationsWaterResidualsByPressure: saturationsWaterResidualsBySaturationsWater;

    Real sumOfUnnormalizedProbabilities = 0;
    if (currentState.isAPressure) {
        // for pressures there's the probability to stay at the same place
        const CellIndex meToMyself = pressureToTransmissibilityIndex(currentState.cell, currentState.cell, numberOfRows);

        const Real correspondingEntryOfAMatrix = 1 - meToMyself(pressureResidualsByPressures);
        const Real correspondingEntryOfBVector = b(currentState.cell, iAmAPressure);
        const Real unnormalizedProbabilityOfStayingHere = std::abs(correspondingEntryOfAMatrix);

        if (unnormalizedProbabilityOfStayingHere > 0) {
            Candidate candidate = {currentState.cell, iAmAPressure, correspondingEntryOfAMatrix,
                                   correspondingEntryOfBVector};
            candidates.push_back(std::move(candidate));
            candidateUnnormalizedProbabilities.push_back(unnormalizedProbabilityOfStayingHere);
        }

        sumOfUnnormalizedProbabilities += unnormalizedProbabilityOfStayingHere;
    }

    // add neighbors
    for (const auto neighbor: currentState.cell.neighbors(numberOfRows, numberOfCols)) {
        ASSERT(neighbor.i >= 0);
        ASSERT(neighbor.j >= 0);
        ASSERT(neighbor.i < numberOfRows);
        ASSERT(neighbor.j < numberOfCols);

        const CellIndex neighborToMe = pressureToTransmissibilityIndex(neighbor, currentState.cell, numberOfRows);

        const Real correspondingEntryOfAMatrixPressure = -neighborToMe(pressureResidualsDerived);
        const Real correspondingEntryOfBVectorPressure = b(neighbor, neighborIsAPressure);
        const Real candidateUnnormalizedProbabilityPressure = std::abs(correspondingEntryOfAMatrixPressure);

        const Real correspondingEntryOfAMatrixSaturationWater = -neighborToMe(saturationWaterResidualsDerived);
        const Real correspondingEntryOfBVectorSaturationWater = b(neighbor, neighborIsASaturation);
        const Real candidateUnnormalizedProbabilitySaturationWater = std::abs(correspondingEntryOfAMatrixSaturationWater);

        if (candidateUnnormalizedProbabilityPressure > 0) {
            Candidate candidatePressure = {neighbor, iAmAPressure, correspondingEntryOfAMatrixPressure,
                                           correspondingEntryOfBVectorPressure};
            candidates.push_back(std::move(candidatePressure));
            candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilityPressure);
        }


        if (candidateUnnormalizedProbabilitySaturationWater > 0) {
            Candidate candidateSaturationWater = {neighbor, iAmASaturation, correspondingEntryOfAMatrixSaturationWater,
                                                  correspondingEntryOfBVectorSaturationWater};
            candidates.push_back(std::move(candidateSaturationWater));
            candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilitySaturationWater);
        }

        sumOfUnnormalizedProbabilities += candidateUnnormalizedProbabilityPressure + candidateUnnormalizedProbabilitySaturationWater;

    }

    const Candidate absorptionCandidate = {CellIndex::invalidCell(), false, NAN, NAN };
    const Real unnormalizedAbsorptionProbability = std::abs(b(currentState.cell, currentState.isAPressure));

    if (unnormalizedAbsorptionProbability > 0) {
        candidates.push_back(absorptionCandidate);
        candidateUnnormalizedProbabilities.push_back(unnormalizedAbsorptionProbability);
    }

    for (const auto& candidate: candidates) {
        ASSERT(candidate.cellIndex.i >= -1);
        ASSERT(candidate.cellIndex.j >= -1);
        ASSERT(candidate.cellIndex.i < numberOfRows);
        ASSERT(candidate.cellIndex.j < numberOfCols);
    }

    std::discrete_distribution<int> choose(candidateUnnormalizedProbabilities.cbegin(), candidateUnnormalizedProbabilities.cend());
    const int nextIndex = choose(rng);

    if (candidates.size() == 0) {
        stillInTheSameTimestep = false;
        return stillInTheSameTimestep;
    }
    const Candidate chosenCandidate = candidates.at(nextIndex);
    ASSERT(candidates.size() == candidateUnnormalizedProbabilities.size());

    const bool willBeAbsorbed = chosenCandidate.cellIndex == absorptionCandidate.cellIndex;

    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;

    if (willBeAbsorbed) {
        stillInTheSameTimestep = false;
        return stillInTheSameTimestep;
    }

    stillInTheSameTimestep = true;
    if (!(currentState.isAPressure && chosenCandidate.isAPressure)) {
        ++currentState.currentTimelevel;
        stillInTheSameTimestep = false;
    }
    // update W (pg. 6199, top)
    currentState.W *= sumOfUnnormalizedProbabilities * chosenCandidate.correspondingEntryOfAMatrix;
    // update D (pg. 6199, top)
    currentState.D += currentState.W * chosenCandidate.correspondingEntryOfBVector;




    ASSERT(currentState.cell.i < numberOfRows);
    ASSERT(currentState.cell.j < numberOfCols);
    ASSERT(currentState.cell.i >= 0);
    ASSERT(currentState.cell.j >= 0);

    return stillInTheSameTimestep;
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
            initialState.D = initialState.W * b(neighborOrMyself, wantAPressure);
            initialState.parameterIndex = parameterIndex;

            randomWalks.push_back(std::move(initialState));
        }
    }

    return randomWalks;



}
