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
                     const SparseMatrix& saturationsWaterResidualsBySaturationsWater, const Real convergenceFactor,
                     const int numberOfRows, const int numberOfCols, Rng& rng) {

    struct Candidate {
        const CellIndex cellIndex;
        const bool isAPressure;
        const Real correspondingEntryOfAMatrix;
        const Real correspondingEntryOfBVector;
    };



    std::vector<Real> candidateUnnormalizedProbabilities;
    std::vector<Candidate> candidates;




    const auto cleanupAndReturn = [&] (const bool stillInTheSameTimestep) {
        return stillInTheSameTimestep;
    };


    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(currentState.D));

    bool stillInTheSameTimestep;
    if (currentState.cell == CellIndex::invalidCell()) {
        stillInTheSameTimestep = false;
        return cleanupAndReturn(stillInTheSameTimestep);
    }

    ASSERT(currentState.cell.i < numberOfRows);
    ASSERT(currentState.cell.j < numberOfCols);
    ASSERT(currentState.cell.i >= 0);
    ASSERT(currentState.cell.j >= 0);







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

        const Real correspondingEntryOfAMatrix = (1 - meToMyself(pressureResidualsByPressures)) * convergenceFactor;
        const Real correspondingEntryOfBVector = b(currentState.cell, iAmAPressure) * convergenceFactor;
        const Real unnormalizedProbabilityOfStayingHere = std::abs(1 - meToMyself(pressureResidualsByPressures));

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

        const Real correspondingEntryOfAMatrixPressure = -neighborToMe(pressureResidualsDerived) * convergenceFactor;
        ASSERT(std::isfinite(correspondingEntryOfAMatrixPressure));
        const Real correspondingEntryOfBVectorPressure = b(neighbor, neighborIsAPressure) * convergenceFactor;
        ASSERT(std::isfinite(correspondingEntryOfBVectorPressure));

        const Real candidateUnnormalizedProbabilityPressure = std::abs(neighborToMe(pressureResidualsDerived));

        const Real correspondingEntryOfAMatrixSaturationWater = -neighborToMe(saturationWaterResidualsDerived) * convergenceFactor;
        ASSERT(std::isfinite(correspondingEntryOfAMatrixSaturationWater));
        const Real correspondingEntryOfBVectorSaturationWater = b(neighbor, neighborIsASaturation) * convergenceFactor;
        ASSERT(std::isfinite(correspondingEntryOfBVectorSaturationWater));
        const Real candidateUnnormalizedProbabilitySaturationWater = std::abs(neighborToMe(saturationWaterResidualsDerived));

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

    constexpr bool alwaysAbsorbAtDrill = false;
    const Real unnormalizedAbsorptionProbability = alwaysAbsorbAtDrill && currentState.isAPressure && currentState.cell == findDrillCell(numberOfRows, numberOfCols)? 1e8: std::abs(b(currentState.cell, currentState.isAPressure));
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


    if (currentState.isAPressure) {
        stillInTheSameTimestep = true;
    } else {
        ++currentState.currentTimelevel;
        stillInTheSameTimestep = false;
    }

    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;

    const bool willBeAbsorbed = chosenCandidate.cellIndex == absorptionCandidate.cellIndex;
    if (willBeAbsorbed) {
        stillInTheSameTimestep = false;
        return cleanupAndReturn(stillInTheSameTimestep);
    }


    // update W (pg. 6199, top)
    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(sumOfUnnormalizedProbabilities));
    ASSERT(std::isfinite(chosenCandidate.correspondingEntryOfAMatrix));
    currentState.W *= sumOfUnnormalizedProbabilities * convergenceFactor * sign(chosenCandidate.correspondingEntryOfAMatrix);

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


    return cleanupAndReturn(stillInTheSameTimestep);
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


void addNewRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters,
                       const int currentTimelevel, const BVectorSurrogate& b, const CMatrixSurrogate& c,
                       std::vector<RandomWalkState>& randomWalks, Rng& rng) {
    constexpr bool justOneRandomWalk = false;
    if (justOneRandomWalk) {
        RandomWalkState justBeginning;
        justBeginning.isAPressure = true;
        justBeginning.cell = {0, 0};
        justBeginning.currentTimelevel = 0;
        justBeginning.W = c(justBeginning.cell, justBeginning.cell, justBeginning.isAPressure);
        justBeginning.D = justBeginning.W * b(justBeginning.cell, justBeginning.isAPressure);
        justBeginning.parameterIndex = 0;
        randomWalks.push_back(justBeginning);
        return;
    }


    constexpr bool initializeJustAtBeginning = false;

    if (initializeJustAtBeginning) {
        LOGGER->info("Initializing random walks");
    } else {
        LOGGER->info("Adding new random walks");
    }

    if (initializeJustAtBeginning && currentTimelevel > 0) {
        return;
    }
    const int numberOfRandomWalksToAdd = 1;



    for (int parameterIndex = 0; parameterIndex <  numberOfParameters; ++parameterIndex) {
        const CellIndex cell = CellIndex::fromLinearIndex(parameterIndex, numberOfRows);

        std::vector<RandomWalkState> candidates;
        std::vector<Real> probabilities;

        std::vector<CellIndex> neighborsAndMyself;

        neighborsAndMyself.push_back(cell);

        constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
              CellIndex::Direction::EAST, CellIndex::Direction::WEST, CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH
        };

        for (const auto direction: directionsToCheck) {
            if (cell.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                neighborsAndMyself.push_back(cell.neighbor(direction));
            }
        }

        constexpr bool justPressure = false;
        constexpr std::array<bool, 2 - justPressure> pressureRequired = {true};

        const Real cNorm = sumOfAbsEntries(c.pressureResidualsByLogPermeability.col(parameterIndex)) + sumOfAbsEntries(c.saturationWaterResidualsByLogPermeability.col(parameterIndex));

        for (const bool wantAPressure: pressureRequired) {

            for (const CellIndex& neighborOrMyself: neighborsAndMyself) {

                const Real cValue = c(neighborOrMyself, cell, wantAPressure);


                const Real prob = std::abs(cValue) / cNorm;
                if (prob < 1e-12) {
                    continue;
                }

                RandomWalkState initialState;
                initialState.cell = neighborOrMyself;
                initialState.isAPressure = wantAPressure;
                initialState.currentTimelevel = currentTimelevel;
                initialState.W =  cNorm * sign(cValue);

                constexpr bool outputInitialization = false;
                if (outputInitialization) {
                    LOGGER->debug("Initialization: want pressure = {}", wantAPressure);
                    LOGGER->debug("Initialization: neighborOrMyself = {}", neighborOrMyself);
                    LOGGER->debug("Initialization: cell = {}", cell);
                    LOGGER->debug("Initialization c-value = {}", c(neighborOrMyself, cell, wantAPressure));
                }

                initialState.D = initialState.W * b(neighborOrMyself, wantAPressure);
                initialState.parameterIndex = parameterIndex;

                candidates.push_back(std::move(initialState));
                probabilities.push_back(prob);
            }
        }

        std::discrete_distribution<int> choose(probabilities.cbegin(), probabilities.cend());

        for (int i = 0; i < numberOfRandomWalksToAdd; ++i) {
            randomWalks.push_back(candidates[choose(rng)]);
        }
    }

    if (initializeJustAtBeginning) {
        LOGGER->info("Finished initializing {} random walks", randomWalks.size());
    } else {
        LOGGER->info("Finished adding {} random walks", numberOfRandomWalksToAdd);
    }



}
