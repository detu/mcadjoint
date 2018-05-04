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


    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(currentState.D));

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
        ASSERT(std::isfinite(correspondingEntryOfAMatrixPressure));
        const Real correspondingEntryOfBVectorPressure = b(neighbor, neighborIsAPressure);
        ASSERT(std::isfinite(correspondingEntryOfBVectorPressure));

        const Real candidateUnnormalizedProbabilityPressure = std::abs(correspondingEntryOfAMatrixPressure);

        const Real correspondingEntryOfAMatrixSaturationWater = -neighborToMe(saturationWaterResidualsDerived);
        ASSERT(std::isfinite(correspondingEntryOfAMatrixSaturationWater));
        const Real correspondingEntryOfBVectorSaturationWater = b(neighbor, neighborIsASaturation);
        ASSERT(std::isfinite(correspondingEntryOfBVectorSaturationWater));
        const Real candidateUnnormalizedProbabilitySaturationWater = std::abs(correspondingEntryOfAMatrixSaturationWater);

        if (candidateUnnormalizedProbabilityPressure > 0) {
            Candidate candidatePressure = {neighbor, iAmAPressure, correspondingEntryOfAMatrixPressure,
                                           correspondingEntryOfBVectorPressure};
            candidates.push_back(std::move(candidatePressure));
            candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilityPressure);
        }

        ASSERT(std::isfinite(candidateUnnormalizedProbabilityPressure));


        if (candidateUnnormalizedProbabilitySaturationWater > 0) {
            Candidate candidateSaturationWater = {neighbor, iAmASaturation, correspondingEntryOfAMatrixSaturationWater,
                                                  correspondingEntryOfBVectorSaturationWater};
            candidates.push_back(std::move(candidateSaturationWater));
            candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilitySaturationWater);
        }

        ASSERT(std::isfinite(candidateUnnormalizedProbabilitySaturationWater));

        sumOfUnnormalizedProbabilities += candidateUnnormalizedProbabilityPressure + candidateUnnormalizedProbabilitySaturationWater;

    }

    const Candidate absorptionCandidate = {CellIndex::invalidCell(), false, 0, 0 };
    const Real unnormalizedAbsorptionProbability = std::abs(b(currentState.cell, currentState.isAPressure));
    ASSERT(unnormalizedAbsorptionProbability > 0 == (currentState.cell == findDrillCell(numberOfRows, numberOfCols) && currentState.isAPressure));
    ASSERT(std::isfinite(unnormalizedAbsorptionProbability));
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

    ASSERT(candidates.size() > 0);
    const Candidate chosenCandidate = candidates.at(nextIndex);
    ASSERT(candidates.size() == candidateUnnormalizedProbabilities.size());


    if (!currentState.isAPressure) {
        ++currentState.currentTimelevel;
        stillInTheSameTimestep = false;
    }

    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;

    const bool willBeAbsorbed = chosenCandidate.cellIndex == absorptionCandidate.cellIndex;
    if (willBeAbsorbed) {
        stillInTheSameTimestep = false;
        return stillInTheSameTimestep;
    }

    stillInTheSameTimestep = true;

    // update W (pg. 6199, top)
    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(sumOfUnnormalizedProbabilities));
    ASSERT(std::isfinite(chosenCandidate.correspondingEntryOfAMatrix));
    currentState.W *= sumOfUnnormalizedProbabilities * chosenCandidate.correspondingEntryOfAMatrix;
    constexpr bool outputW = false;
    if (outputW) {
        LOGGER->debug("Will be pressure = {}", chosenCandidate.isAPressure);
        LOGGER->debug("Will be cell index = {} ",  chosenCandidate.cellIndex);
        LOGGER->debug("W = {}", currentState.W);
        LOGGER->debug("Sum of unnormalized Probs = {}", sumOfUnnormalizedProbabilities);
        LOGGER->debug("A entry = {}", chosenCandidate.correspondingEntryOfAMatrix);
    }
    ASSERT(std::isfinite(currentState.W));
    // update D (pg. 6199, top)
    ASSERT(std::isfinite(currentState.D));
    ASSERT(std::isfinite(chosenCandidate.correspondingEntryOfBVector));
    currentState.D += currentState.W * chosenCandidate.correspondingEntryOfBVector;
    ASSERT(std::isfinite(currentState.D));





    ASSERT(currentState.cell.i < numberOfRows);
    ASSERT(currentState.cell.j < numberOfCols);
    ASSERT(currentState.cell.i >= 0);
    ASSERT(currentState.cell.j >= 0);

    return stillInTheSameTimestep;
}

void logStatisticsAboutRandomWalks(const std::vector<RandomWalkState>& randomWalks) {
    LOGGER->info("We have {} random walks", randomWalks.size());

    int absorbedWalks = 0;
    int walksInPressure = 0;
    int walksInSaturation = 0;

    const CellIndex invalidCell = CellIndex::invalidCell();
    for (const RandomWalkState& randomWalk: randomWalks) {
        if (randomWalk.cell == invalidCell) {
            ++absorbedWalks;
        } else if (randomWalk.isAPressure) {
            ++walksInPressure;
        } else {
            ++walksInSaturation;
        }
    }

    LOGGER->info("{} walks are currently in pressure cells", walksInPressure);
    LOGGER->info("{} walks are currently in saturation cells", walksInSaturation);
    LOGGER->info("{} walks are absorbed", absorbedWalks);

    ASSERT(absorbedWalks + walksInPressure + walksInSaturation == randomWalks.size());

}


std::vector<RandomWalkState> initializeRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters, const BVectorSurrogate& b, const CMatrixSurrogate& c) {
    constexpr bool justOneRandomWalk = false;
    if (justOneRandomWalk) {
        RandomWalkState justBeginning;
        justBeginning.isAPressure = true;
        justBeginning.cell = {0, 0};
        justBeginning.currentTimelevel = 0;
        justBeginning.W = c(justBeginning.cell, justBeginning.cell, justBeginning.isAPressure);
        justBeginning.D = justBeginning.W * b(justBeginning.cell, justBeginning.isAPressure);
        justBeginning.parameterIndex = 0;
        return {justBeginning};
    }
    const int numberOfRandomWalksPerPressureCell = 10;

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
            LOGGER->debug("Initialization: want pressure = {}", wantAPressure);
            LOGGER->debug("Initialization: neighborOrMyself = {}", neighborOrMyself);
            LOGGER->debug("Initialization: cell = {}", cell);
            LOGGER->debug("Initialization c-value = {}", c(neighborOrMyself, cell, wantAPressure));
            initialState.D = initialState.W * b(neighborOrMyself, wantAPressure);
            initialState.parameterIndex = parameterIndex;

            if (initialState.W != 0) {
                randomWalks.push_back(std::move(initialState));
            }
        }
    }

    return randomWalks;



}
