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
#include "pressure.hpp"
#include "utils.hpp"


static inline int cellIndexToBIndex(const CellIndex& cellIndex, const bool isAPressure, const int numberOfRows, const int numberOfCols) {
    int linearIndex = cellIndex.linearIndex(numberOfRows);

    if (!isAPressure) {
        linearIndex += numberOfCols * numberOfRows;
    }

    return linearIndex;
}

bool transitionState(RandomWalkState& currentState, ConstVectorRef b,
                     const SparseMatrix& pressureResidualsByPressures,
                     const SparseMatrix& pressureResidualsBySaturationsWater,
                     const SparseMatrix& saturationsWaterResidualsByPressure,
                     const SparseMatrix& saturationsWaterResidualsBySaturationsWater, const int numberOfRows,
                     const int numberOfCols, Rng& rng) {

    struct Candidate {
        const CellIndex cellIndex;
        const bool isAPressure;
        const Real correspondingEntryOfAMatrix;
        const Real correspondingEntryOfBVector;
    };


    constexpr Candidate absorptionCandidate = {CellIndex::invalidCell(), false, 0, 0};


    constexpr bool surelyNotStillInTheSameTimestep  = false;

    if (currentState.cell == absorptionCandidate.cellIndex) {
        return surelyNotStillInTheSameTimestep;
    }

    bool stillInTheSameTimestep;


    std::vector<Real> candidateUnnormalizedProbabilities;
    std::vector<Candidate> candidates;


    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(currentState.D));



    ASSERT(currentState.cell.i < numberOfRows);
    ASSERT(currentState.cell.j < numberOfCols);
    ASSERT(currentState.cell.i >= 0);
    ASSERT(currentState.cell.j >= 0);


    const SparseMatrix& pressureResidualsDerived = (currentState.isAPressure? pressureResidualsByPressures: pressureResidualsBySaturationsWater);
    const SparseMatrix& saturationWaterResidualsDerived = currentState.isAPressure? saturationsWaterResidualsByPressure: saturationsWaterResidualsBySaturationsWater;


    Real sumOfUnnormalizedProbabilities = 0;


    for (const CellIndex& target: currentState.cell.neighborsAndMyself(numberOfRows, numberOfCols)) {
        ASSERT(target.i >= 0);
        ASSERT(target.j >= 0);
        ASSERT(target.i < numberOfRows);
        ASSERT(target.j < numberOfCols);


        for (const bool targetIsPressure: std::array<bool, 2>({true, false})) {
            const SparseMatrix& correspondingDerivative = targetIsPressure? pressureResidualsDerived
                                                                            : saturationWaterResidualsDerived;


            const CellIndex neighborToMe = pressureToTransmissibilityIndex(target, currentState.cell, numberOfRows);

            const Real correspondingEntryOfAMatrix = Real(targetIsPressure && currentState.isAPressure && target == currentState.cell) - neighborToMe(correspondingDerivative);
            //log()->debug("from ({}, {}) to ({}, {})", currentState.cell, currentState.isAPressure, target, targetIsPressure);
            //log()->debug("corresponding entry of a matrix = {}", correspondingEntryOfAMatrix);
            ASSERT(std::isfinite(correspondingEntryOfAMatrix));
            const Real correspondingEntryOfBVector = b(cellIndexToBIndex(target, targetIsPressure, numberOfRows, numberOfCols));
            ASSERT(std::isfinite(correspondingEntryOfBVector));

            const Real candidateUnnormalizedProbability = std::abs(correspondingEntryOfAMatrix);


            if (candidateUnnormalizedProbability > 0) {
                Candidate candidatePressure = {target, targetIsPressure, correspondingEntryOfAMatrix,
                                               correspondingEntryOfBVector};
                candidates.push_back(std::move(candidatePressure));
                candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbability);
            }

            ASSERT(std::isfinite(candidateUnnormalizedProbability));


            sumOfUnnormalizedProbabilities += candidateUnnormalizedProbability;
        }

    }



    {

        constexpr bool alwaysAbsorbAtDrill = false;
        const Real unnormalizedAbsorptionProbability = alwaysAbsorbAtDrill && currentState.isAPressure &&
                                                       currentState.cell == findDrillCell(numberOfRows, numberOfCols)
                                                       ? 1e8: std::abs(b(cellIndexToBIndex(currentState.cell, currentState.isAPressure, numberOfRows, numberOfCols)));
        ASSERT(std::isfinite(unnormalizedAbsorptionProbability));

        if (unnormalizedAbsorptionProbability > 0 && false) {
            candidates.push_back(absorptionCandidate);
            candidateUnnormalizedProbabilities.push_back(unnormalizedAbsorptionProbability);
        }
    }


    for (const auto& candidate: candidates) {
        ASSERT(candidate.cellIndex.i >= -1);
        ASSERT(candidate.cellIndex.j >= -1);
        ASSERT(candidate.cellIndex.i < numberOfRows);
        ASSERT(candidate.cellIndex.j < numberOfCols);
    }


    //ASSERT(!candidates.empty());
    if (candidates.empty()) {
        ASSERT(candidateUnnormalizedProbabilities.empty());
        candidates.push_back(absorptionCandidate);
        candidateUnnormalizedProbabilities.push_back(1);
    }

    std::discrete_distribution<int> choose(candidateUnnormalizedProbabilities.cbegin(), candidateUnnormalizedProbabilities.cend());
    const int nextIndex = choose(rng);

    const Candidate chosenCandidate = candidates[nextIndex];
    ASSERT(candidates.size() == candidateUnnormalizedProbabilities.size());


    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;


    const bool willBeAbsorbed = chosenCandidate.cellIndex == absorptionCandidate.cellIndex;
    if (willBeAbsorbed) {
        currentState.W = NAN;
        return surelyNotStillInTheSameTimestep;
    }


    if (currentState.isAPressure) {
        stillInTheSameTimestep = true;
    } else {
        ++currentState.currentTimelevel;
        stillInTheSameTimestep = false;
    }






    // update W (pg. 6199, top)
    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(sumOfUnnormalizedProbabilities));
    ASSERT(std::isfinite(chosenCandidate.correspondingEntryOfAMatrix));

    const Real probabilityOfChoosingThisCandidate = candidateUnnormalizedProbabilities[nextIndex] / sumOfUnnormalizedProbabilities;
    ASSERT(sumOfUnnormalizedProbabilities > 0);
    ASSERT(probabilityOfChoosingThisCandidate > 0);
    currentState.W *= chosenCandidate.correspondingEntryOfAMatrix / probabilityOfChoosingThisCandidate;
    ASSERT(std::isfinite(currentState.W));


    /*constexpr bool outputW = false;
    if (outputW) {
        log()->debug("Will be pressure = {}", chosenCandidate.isAPressure);
        log()->debug("Will be cell index = {} ",  chosenCandidate.cellIndex);
        log()->debug("W = {}", currentState.W);
        log()->debug("Sum of unnormalized Probs = {}", sumOfUnnormalizedProbabilities);
        log()->debug("A entry = {}", chosenCandidate.correspondingEntryOfAMatrix);
    }*/

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
    log()->info("We have {} random walks", randomWalks.size());

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

    log()->info("{} walks are currently in pressure cells", walksInPressure);
    log()->info("{} walks are currently in saturation cells", walksInSaturation);
    log()->info("{} walks are absorbed", absorbedWalks);

    ASSERT(absorbedWalks + walksInPressure + walksInSaturation == randomWalks.size());

}

static inline CellIndex cellIndexToCIndex(const CellIndex& cellIndex, const CellIndex& parameterIndex, const bool isAPressure, const int numberOfRows, const int numberOfCols) {
    CellIndex derivativeIndex = pressureToTransmissibilityIndex(cellIndex, parameterIndex, numberOfRows);
    if (!isAPressure) {
        derivativeIndex.i += numberOfCols * numberOfRows;
    }

    return derivativeIndex;
}
void addNewRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters,
                       const int currentTimelevel, ConstVectorRef b, const SparseMatrix& c,
                       std::vector<RandomWalkState>& randomWalks, Rng& rng) {

    constexpr bool initializeJustAtBeginning = false;

    if (initializeJustAtBeginning) {
        log()->info("Initializing random walks");
    } else {
        log()->info("Adding new random walks");
    }

    if (initializeJustAtBeginning && currentTimelevel > 0) {
        return;
    }
    const int numberOfRandomWalksToAdd = 1000;



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

        constexpr std::array<bool, 2> pressureRequired = {true, false};

        const Real cNorm = sumOfAbsEntries(c.col(parameterIndex));

        if (cNorm < 1e-12) {
            continue;
        }
        for (const bool wantAPressure: pressureRequired) {

            for (const CellIndex& neighborOrMyself: neighborsAndMyself) {

                const Real cValue = cellIndexToCIndex(neighborOrMyself, cell, wantAPressure, numberOfRows, numberOfCols)(c);


                const Real prob = std::abs(cValue) / cNorm;
                if (prob < 1e-12) {
                    continue;
                }

                RandomWalkState initialState;
                initialState.cell = neighborOrMyself;
                initialState.isAPressure = wantAPressure;
                initialState.currentTimelevel = currentTimelevel;
                initialState.W =  cValue / prob;

                initialState.D = initialState.W * b(cellIndexToBIndex(neighborOrMyself, wantAPressure, numberOfRows, numberOfCols));
                initialState.parameterIndex = parameterIndex;

                candidates.push_back(std::move(initialState));
                probabilities.push_back(prob);
            }
        }

        for (int candidateIndex = 0; candidateIndex < int(candidates.size()); ++candidateIndex) {
            ASSERT(probabilities[candidateIndex] > 0 && probabilities[candidateIndex] <= 1);
            const int copiesToAddOfThisCandidate = std::round(numberOfRandomWalksToAdd * probabilities[candidateIndex]);
            const RandomWalkState& candidate = candidates[candidateIndex];

            for (int i = 0; i < copiesToAddOfThisCandidate; ++i) {
                randomWalks.push_back(candidate);
            }

        }

    }

    if (initializeJustAtBeginning) {
        log()->info("Finished initializing {} random walks", randomWalks.size());
    } else {
        log()->info("Finished adding {} random walks", numberOfRandomWalksToAdd);
    }



}

