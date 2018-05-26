//
// Created by stefanow on 4/27/18.
//

#include "sensitivity.hpp"
#include <stefCommonHeaders/dev.hpp>
#include "specialCells.hpp"
#include "simulationState.hpp"
#include "forward.hpp"
#include "pressure.hpp"
#include "derivativesForAdjoint.hpp"
#include "adjoint.hpp"
#include "logging.hpp"
#include "regularization.hpp"
#include "dumpToMatFile.hpp"
#include "utils.hpp"
#include "adjointOptions.hpp"

#include <cmath>
#include <list>

Real computeContributionToCost(const FixedParameters& parameters, const SimulationState& currentSimulationState) {

    const int numberOfRows = currentSimulationState.saturationsWater.rows();
    const int numberOfCols = currentSimulationState.saturationsWater.cols();
    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);

    const Real computedPressureAtDrill = drillCell(currentSimulationState.pressures.map);
    log()->info("computedPressureAtDrill = {}", computedPressureAtDrill);
    const Real measuredPressureAtDrill = parameters.overPressureDrill(currentSimulationState.time);
    log()->info("measured pressure at drill = {}", measuredPressureAtDrill);
    return std::pow(measuredPressureAtDrill - computedPressureAtDrill, 2);
}

SensitivityAndCost computeSensitivityAndCostTraditional(const FixedParameters& params, ConstMatrixRef permeabilities,
                          ConstMatrixRef logPermeabilities) {
    log()->info("Computing sensitivity by directly solving the traditional system");
    const int numberOfRows = params.initialPermeabilities.rows();
    const int numberOfCols = params.initialPermeabilities.cols();
    const int numberOfParameters = numberOfRows * numberOfCols;

    const int numberOfTimesteps = params.maxNumberOfTimesteps;

    SimulationState simulationState(numberOfRows, numberOfCols);
    const int stateSize = 2 * numberOfCols * numberOfRows;
    Matrix adjointMatrix = Matrix::Zero(stateSize * numberOfTimesteps, stateSize * numberOfTimesteps);
    Vector adjointRhs = Vector::Zero(stateSize * numberOfTimesteps);
    Matrix completeC = Matrix::Zero(stateSize * numberOfTimesteps, numberOfParameters);

    bool brokeThrough = false;
    Real cost = 0;
    for (int currentTimelevel = 0; currentTimelevel < numberOfTimesteps && !brokeThrough; ++currentTimelevel) {
        brokeThrough = stepForwardAndAdjointProblemTraditional(params, permeabilities, currentTimelevel, simulationState, adjointMatrix, adjointRhs, completeC);

        cost += computeContributionToCost(params, simulationState);
        log()->info("time = {}", simulationState.time);
        //if (!params.traditionalMinimization) {
            dumpThis("adjointMatrixTrad", adjointMatrix);
            dumpThis("adjointRhsTrad", adjointRhs);
            dumpThis("completeCTrad", completeC);
        //}
        ASSERT(allFinite(adjointMatrix));
        ASSERT(allFinite(adjointRhs));
    }
    log()->debug("broke through? {}", brokeThrough);

    const Vector adjoint = adjointMatrix.householderQr().solve(adjointRhs);
    //if (!params.traditionalMinimization) {
        dumpThis("adjointTrad", adjoint);
    //}

    const Vector sensitivities = -completeC.transpose() * adjoint;


    SensitivityAndCost sensitivityAndCost;
    sensitivityAndCost.sensitivity.resizeLike(sensitivities);

    sensitivityAndCost.cost = cost;
    sensitivityAndCost.sensitivity = sensitivities;
    applyRegularizationIfEnabled(params, logPermeabilities, sensitivityAndCost);

    dumpThis("sensitivityTrad", sensitivityAndCost.sensitivity);

    return sensitivityAndCost;


}

static inline CellIndex getSymmetricCounterpart(const CellIndex& cellIndex, const int numberOfRows, const int numberOfCols) {
    return {numberOfCols - 1 - cellIndex.j, numberOfRows - 1 - cellIndex.i};
}

static inline int getSymmetricCounterpart(const int linearIndex, const int numberOfRows, const int numberOfCols) {
    return getSymmetricCounterpart(CellIndex::fromLinearIndex(linearIndex, numberOfRows), numberOfRows, numberOfCols).linearIndex(numberOfRows);
}


SensitivityAndCost computeSensitivityAndCost(const FixedParameters& params, ConstMatrixRef permeabilities,
                                             ConstMatrixRef logPermeabilities, Rng& rng) {
    log()->info("Computing sensitivity using MC");

    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();

    SimulationState simulationState(numberOfRows, numberOfCols);
    std::list<RandomWalkState> randomWalks;
    std::list<RandomWalkState> antitheticRandomWalks;


    Eigen::VectorXi numberOfRandomWalksPerParameter = Eigen::VectorXi::Zero(numberOfParameters);
    Eigen::VectorXi numberOfRandomWalksPerParameterAntithetic = Eigen::VectorXi::Zero(numberOfParameters);
    PreciseVector sensitivities = PreciseVector::Zero(numberOfParameters);
    PreciseVector sensitivitiesAntithetic = PreciseVector::Zero(numberOfParameters);

    bool breakthroughHappened = false;
    int currentTimeLevel = 0;

    Real cost = 0;

    do {
        log()->info("-----------------------------------");
        log()->info("time = {}", simulationState.time);
        breakthroughHappened = stepForwardAndAdjointProblem(params, permeabilities, currentTimeLevel, params.numberOfRandomWalksToAdd, simulationState,
                                                            randomWalks, antitheticRandomWalks, rng);
        removeAbsorbedStates(randomWalks, numberOfRandomWalksPerParameter, sensitivities);
        if (params.enableAntitheticSampling) {
            log()->info("Antithetic sampling enabled");
            removeAbsorbedStates(antitheticRandomWalks, numberOfRandomWalksPerParameterAntithetic, sensitivitiesAntithetic);
        } else {
            log()->info("Antithetic sampling disabled");
        }


        const Real contributionToCost = computeContributionToCost(params, simulationState);
        cost += contributionToCost;
        ++currentTimeLevel;
    } while (!breakthroughHappened && simulationState.time < params.finalTime && currentTimeLevel < params.maxNumberOfTimesteps);



    using RandomWalkIterator = typename std::list<RandomWalkState>::iterator;
    RandomWalkIterator randomWalkIterator = randomWalks.begin();
    RandomWalkIterator antitheticRandomWalkIterator = antitheticRandomWalks.begin();
    const RandomWalkIterator randomWalkEnd = randomWalks.end();
    const RandomWalkIterator antitheticRandomWalkEnd = antitheticRandomWalks.end();


    while (randomWalkIterator != randomWalkEnd) {
        sensitivities(randomWalkIterator->parameterIndex) += randomWalkIterator->D;
        ++numberOfRandomWalksPerParameter(randomWalkIterator->parameterIndex);
        ++randomWalkIterator;
    }

    if (params.enableAntitheticSampling) {
        while (antitheticRandomWalkIterator != antitheticRandomWalkEnd) {
            sensitivitiesAntithetic(antitheticRandomWalkIterator->parameterIndex) += antitheticRandomWalkIterator->D;
            ++numberOfRandomWalksPerParameterAntithetic(antitheticRandomWalkIterator->parameterIndex);
            ++antitheticRandomWalkIterator;
        }
    }

    dumpThis("numberOfRandomWalksPerParameter", numberOfRandomWalksPerParameter.cast<Real>());
    dumpThis("numberOfRandomWalksPerParameterAntithetic", numberOfRandomWalksPerParameterAntithetic.cast<Real>());



    sensitivities.array() /= numberOfRandomWalksPerParameter.array().cwiseMax(1).cast<WiderReal>();
    /*if (!initializeJustAtBeginning) {
        sensitivities.array() *= currentTimeLevel;
    }*/

    if (params.enableAntitheticSampling) {
        sensitivitiesAntithetic.array() /= numberOfRandomWalksPerParameterAntithetic.array().cwiseMax(1).cast<WiderReal>();
        /*if (!initializeJustAtBeginning) {
            sensitivitiesAntithetic.array() *= currentTimeLevel;
        }*/
    }



    Real antitheticPart = 0;
    if (params.enableAntitheticSampling) {
        antitheticPart = 0.5;
    }

    SensitivityAndCost sensitivityAndCost;
    sensitivityAndCost.sensitivity.resizeLike(sensitivities);

    sensitivityAndCost.cost = cost;

    log()->info("sensitivities norm = {}", sensitivities.norm());


    if (params.enableAntitheticSampling) {
        sensitivityAndCost.sensitivity =
              ((1.0 - antitheticPart) * sensitivities + antitheticPart * sensitivitiesAntithetic).cast<Real>();

        log()->info("antithetic sensitivities norm = {}", sensitivitiesAntithetic.norm());
    } else {
        sensitivityAndCost.sensitivity = sensitivities.cast<Real>();
    }

    sensitivityAndCost.sensitivity *= -1;

    if (params.symmetrizeGradient) {
        log()->info("Symmetrizing gradient");

        Vector copyOfSensitivities = sensitivityAndCost.sensitivity;
        std::nth_element(copyOfSensitivities.data(), copyOfSensitivities.data() + copyOfSensitivities.size()/2, copyOfSensitivities.data() + copyOfSensitivities.size());
        const Real median = copyOfSensitivities(copyOfSensitivities.size() / 2);
        const Real mean = sensitivityAndCost.sensitivity.mean();
        const Real stdDev = std::sqrt(1.0 /(sensitivityAndCost.sensitivity.size() - 1) * (sensitivityAndCost.sensitivity.array() - mean).square().sum());
        for (int linearIndex = 0; linearIndex < sensitivityAndCost.sensitivity.size()/2; ++linearIndex) {
            const Real myValue = sensitivityAndCost.sensitivity(linearIndex);
            const int counterpartIndex = getSymmetricCounterpart(linearIndex, numberOfRows, numberOfCols);
            const Real counterpartValue = sensitivityAndCost.sensitivity(counterpartIndex);



            const Real myDistance = std::abs(myValue - median);
            const Real counterpartDistance = std::abs(counterpartValue - median);

            Real moreRealisticValue = NAN;
            Real smallerDistance = NAN;
            if (myDistance < counterpartDistance) {
                moreRealisticValue = myValue;
                smallerDistance = myDistance;
            } else {
                moreRealisticValue = counterpartValue;
                smallerDistance = counterpartDistance;
            }


            Real smootherValue = NAN;

            if (smallerDistance < stdDev) {
                smootherValue = moreRealisticValue;
            } else {
                smootherValue = 0.5 * median + 0.5 * moreRealisticValue;
            }
            sensitivityAndCost.sensitivity(linearIndex) = smootherValue;
            sensitivityAndCost.sensitivity(counterpartIndex) = smootherValue;
        }
    } else {
        log()->info("Not symmetrizing gradient");
    }

    applyRegularizationIfEnabled(params, logPermeabilities, sensitivityAndCost);
    #ifdef JUST_COMPUTE_ADJOINT
        #pragma message "Just computing adjoint"

    dumpThis("adjointMC", sensitivityAndCost.sensitivity);
        #endif
    dumpThis("sensitivities", sensitivityAndCost.sensitivity);
    log()->debug("Sensitivities =\n{}", sensitivityAndCost.sensitivity);

    return sensitivityAndCost;
}


