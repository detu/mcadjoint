//
// Created by Stefano Weidmann on 17.04.18.
//

#include "adjoint.hpp"
#include "adjointOptions.hpp"
#include "specialCells.hpp"
#include "fast_discrete_distribution.hpp"
#include <stefCommonHeaders/assert.h>
#include <stefCommonHeaders/xoroshiro.h>
#include <iostream>
#include "logging.hpp"
#include "pressure.hpp"
#include "utils.hpp"
#include "pickAnIndex.hpp"
#include <random>


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
                     const int numberOfCols, const Real standardUniformNumber) {


    struct Candidate {
        const CellIndex cellIndex;
        const bool isAPressure;
        const Real correspondingEntryOfAMatrix;
        const Real correspondingEntryOfBVector;
        const bool advancesTime;
    };

    static constexpr Candidate absorptionCandidate = {CellIndex::invalidCell(), false, 0, 0, true};



    if (currentState.cell == absorptionCandidate.cellIndex) {
        return false;
    }

    bool stillInTheSameTimestep;





    std::vector<Real> candidateUnnormalizedProbabilities;
    std::vector<Candidate> candidates;
    WiderReal sumOfUnnormalizedProbabilities = 0;



    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(currentState.D));



    ASSERT(currentState.cell.i < numberOfRows);
    ASSERT(currentState.cell.j < numberOfCols);
    ASSERT(currentState.cell.i >= 0);
    ASSERT(currentState.cell.j >= 0);


    const SparseMatrix& pressureResidualsDerived = (currentState.isAPressure? pressureResidualsByPressures: pressureResidualsBySaturationsWater);
    const SparseMatrix& saturationWaterResidualsDerived = currentState.isAPressure? saturationsWaterResidualsByPressure: saturationsWaterResidualsBySaturationsWater;

    Real correctionFactorForGamma = 1;

    if (correctForGamma) {
        Real colSum = 0;
        const int myLinearIndex = currentState.cell.linearIndex(numberOfRows);
        if (currentState.isAPressure) {
            colSum = saturationWaterResidualsDerived.col(myLinearIndex).cwiseAbs().sum() + pressureResidualsDerived.col(myLinearIndex).cwiseAbs().sum() ;
        } else {
            colSum = pressureResidualsDerived.col(myLinearIndex).cwiseAbs().sum() +
                     saturationWaterResidualsDerived.col(myLinearIndex).cwiseAbs().sum() + 1;
        }

        correctionFactorForGamma = std::max(1.0, colSum);

        const Real correspondingEntryOfA = 1.0 - 1.0 / correctionFactorForGamma;
        const Real correspondingEntryOfB = b(cellIndexToBIndex(currentState.cell, currentState.isAPressure, numberOfRows, numberOfCols)) / correctionFactorForGamma;
        log()->info("Col sum = {}, corrected col sum = {}", colSum, (colSum > 1? 2.0 * (1 - 1.0 / colSum): colSum));
        const bool advancesTime = false;

        const Real candidateUnnormalizedProbability = std::abs(correspondingEntryOfA);
        sumOfUnnormalizedProbabilities += candidateUnnormalizedProbability;

        if (candidateUnnormalizedProbability > 0.0) {
            candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbability);
            candidates.push_back({currentState.cell, currentState.isAPressure, correspondingEntryOfA, correspondingEntryOfB, advancesTime});
        }


    }




    for (const CellIndex& target: currentState.cell.neighborsAndMyself(numberOfRows, numberOfCols)) {
        ASSERT(target.i >= 0);
        ASSERT(target.j >= 0);
        ASSERT(target.i < numberOfRows);
        ASSERT(target.j < numberOfCols);


        for (const bool targetIsPressure: std::array<bool, 2>({true, false})) {
            const SparseMatrix& correspondingDerivative = targetIsPressure? pressureResidualsDerived
                                                                            : saturationWaterResidualsDerived;


            const CellIndex neighborToMe = pressureToTransmissibilityIndex(target, currentState.cell, numberOfRows);

            if (targetIsPressure && currentState.isAPressure && currentState.cell == target) {
                // handled above
                continue;
            }

            const Real correspondingEntryOfAMatrix = -neighborToMe(correspondingDerivative) / correctionFactorForGamma;
            //log()->debug("from ({}, {}) to ({}, {})", currentState.cell, currentState.isAPressure, target, targetIsPressure);
            //log()->debug("corresponding entry of a matrix = {}", correspondingEntryOfAMatrix);
            ASSERT(std::isfinite(correspondingEntryOfAMatrix));
            const Real  correspondingEntryOfBVector =  b(cellIndexToBIndex(target, targetIsPressure, numberOfRows, numberOfCols)) / correctionFactorForGamma;
            ASSERT(std::isfinite(correspondingEntryOfBVector));

            const Real candidateUnnormalizedProbability = std::abs(correspondingEntryOfAMatrix);


            if (candidateUnnormalizedProbability > 0) {
                Candidate candidate = {target, targetIsPressure, correspondingEntryOfAMatrix,
                                               correspondingEntryOfBVector, !currentState.isAPressure};
                candidates.push_back(std::move(candidate));
                candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbability);
            }

            ASSERT(std::isfinite(candidateUnnormalizedProbability));


            sumOfUnnormalizedProbabilities += candidateUnnormalizedProbability;
        }

    }



    {


        const Real unnormalizedAbsorptionProbability = alwaysAbsorbAtDrill && currentState.isAPressure &&
                                                       currentState.cell == findDrillCell(numberOfRows, numberOfCols)
                                                       ? 1e8: std::abs(b(cellIndexToBIndex(currentState.cell, currentState.isAPressure, numberOfRows, numberOfCols)));
        ASSERT(std::isfinite(unnormalizedAbsorptionProbability));

        if (unnormalizedAbsorptionProbability > 0 && enableAbsorption) {
            candidates.push_back(absorptionCandidate);
            candidateUnnormalizedProbabilities.push_back(unnormalizedAbsorptionProbability);
        }
    }

    #ifndef NDEBUG
    for (const auto& candidate: candidates) {
        ASSERT(candidate.cellIndex.i >= -1);
        ASSERT(candidate.cellIndex.j >= -1);
        ASSERT(candidate.cellIndex.i < numberOfRows);
        ASSERT(candidate.cellIndex.j < numberOfCols);
    }
    #endif


    //ASSERT(!candidates.empty());
    if (candidates.empty()) {
        ASSERT(candidateUnnormalizedProbabilities.empty());
        candidates.push_back(absorptionCandidate);
        candidateUnnormalizedProbabilities.push_back(1);
    }



    const int nextIndex = pickAnIndex(candidateUnnormalizedProbabilities, standardUniformNumber);

    const Candidate chosenCandidate = candidates[nextIndex];
    ASSERT(candidates.size() == candidateUnnormalizedProbabilities.size());


    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;


    const bool willBeAbsorbed = chosenCandidate.cellIndex == absorptionCandidate.cellIndex;
    if (willBeAbsorbed) {
        currentState.W = NAN;
        return false;
    }


    if (chosenCandidate.advancesTime) {
        ++currentState.currentTimelevel;
        stillInTheSameTimestep = false;
    } else {
        stillInTheSameTimestep = true;
    }







    // update W (pg. 6199, top)
    ASSERT(std::isfinite(currentState.W));
    ASSERT(std::isfinite(sumOfUnnormalizedProbabilities));
    ASSERT(std::isfinite(chosenCandidate.correspondingEntryOfAMatrix));

    const WiderReal  probabilityOfChoosingThisCandidate = candidateUnnormalizedProbabilities[nextIndex] / sumOfUnnormalizedProbabilities;
    ASSERT(sumOfUnnormalizedProbabilities > 0);
    ASSERT(probabilityOfChoosingThisCandidate > 0);
    currentState.W *= chosenCandidate.correspondingEntryOfAMatrix / probabilityOfChoosingThisCandidate;
    ASSERT(std::isfinite(currentState.W));



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

void logStatisticsAboutRandomWalks(const std::list<RandomWalkState>& randomWalks) {
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

    ASSERT(absorbedWalks + walksInPressure + walksInSaturation == int(randomWalks.size()));

}

static inline CellIndex cellIndexToCIndex(const CellIndex& cellIndex, const CellIndex& parameterIndex, const bool isAPressure, const int numberOfRows, const int numberOfCols) {
    CellIndex derivativeIndex = pressureToTransmissibilityIndex(cellIndex, parameterIndex, numberOfRows);
    if (!isAPressure) {
        derivativeIndex.i += numberOfCols * numberOfRows;
    }

    return derivativeIndex;
}

static inline void addCopiesOfState(const RandomWalkState& toAdd, const Real probability, const int numberOfRandomWalksToAdd, std::list<RandomWalkState>& randomWalksToAddTo) {
    const int numberOfCopiesToAdd = std::round(probability * numberOfRandomWalksToAdd);

    for (int copyIndex = 0; copyIndex < numberOfCopiesToAdd; ++copyIndex) {
        randomWalksToAddTo.push_back(toAdd);
    }
}

void addNewRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters,
                       const int currentTimelevel, const int numberOfRandomWalksToAdd, const bool enableAntitheticSampling, ConstVectorRef b, SparseMatrix c,
                       std::list<RandomWalkState>& randomWalks, std::list<RandomWalkState>& antitheticRandomWalks, Rng& rng) {



    if (initializeJustAtBeginning) {
        log()->info("Initializing random walks");
    } else {
        log()->info("Adding new random walks");
    }

    if (initializeJustAtBeginning && currentTimelevel > 0) {
        return;
    }

    const int numberOfPressures = numberOfCols * numberOfRows;
    if (preferSaturations) {
        log()->info("Prefering saturations with a factor of {}", preferenceForSaturations);
        for (int outerIndex = 0; outerIndex < c.outerSize(); ++outerIndex) {
            for (SparseMatrix::InnerIterator innerIterator(c, outerIndex); bool(innerIterator); ++innerIterator) {
                const bool correspondsToASaturation = innerIterator.row() >= numberOfPressures;

                if (correspondsToASaturation) {
                    innerIterator.valueRef() *= preferenceForSaturations;
                }
            }
        }
    } else {
        log()->info("Not preferring saturations over pressures");
    }

    const Real minimumProbabilityToBeAdded = 0.5 / Real(numberOfRandomWalksToAdd);



    for (int parameterIndex = 0; parameterIndex <  numberOfParameters; ++parameterIndex) {
        const CellIndex cell = CellIndex::fromLinearIndex(parameterIndex, numberOfRows);


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


                const WiderReal prob = std::abs(cValue) / cNorm;
                if (prob < minimumProbabilityToBeAdded) {
                        continue;
                }

                RandomWalkState initialState;
                initialState.cell = neighborOrMyself;
                initialState.isAPressure = wantAPressure;
                initialState.currentTimelevel = currentTimelevel;
                initialState.startingTimeLevel = currentTimelevel;
                initialState.W = ((WiderReal) cValue) / prob;

                initialState.D = initialState.W * b(cellIndexToBIndex(neighborOrMyself, wantAPressure, numberOfRows, numberOfCols));
                initialState.parameterIndex = parameterIndex;

                addCopiesOfState(initialState, prob, numberOfRandomWalksToAdd, randomWalks);

                if (enableAntitheticSampling) {
                    addCopiesOfState(initialState, prob, numberOfRandomWalksToAdd, antitheticRandomWalks);
                }

            }
        }



    }


    if (initializeJustAtBeginning) {
        log()->info("Finished initializing {} random walks", randomWalks.size());
    } else {
        log()->info("Finished adding {} random walks", numberOfRandomWalksToAdd);
    }


}

void removeAbsorbedStates(std::list<RandomWalkState>& randomWalks, Eigen::VectorXi& numberOfRemovedAbsorbedStates, PreciseVector& sumOfDValuesOfAbsorbedStates) {


    log()->info("Clearing absorbed states");
    auto currentState = randomWalks.begin();
    const auto endState = randomWalks.end();

    if (currentState == endState) {
        log()->info("Nothing to clear");

        return;
    }

    constexpr CellIndex cellOfAnAbsorbedState = CellIndex::invalidCell();


    while (currentState != endState) {

        const bool currentStateIsAbsorbed = currentState->cell == cellOfAnAbsorbedState;

        if (currentStateIsAbsorbed) {
            const int parameterIndexOfThisAbsorbedState = currentState->parameterIndex;
            ++numberOfRemovedAbsorbedStates(parameterIndexOfThisAbsorbedState);
            sumOfDValuesOfAbsorbedStates(parameterIndexOfThisAbsorbedState) += currentState->D;
            currentState = randomWalks.erase(currentState);
        } else {
            ++currentState;
        }

    }

    log()->info("Cleared absorbed states");

}