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
#include "antitheticOptions.hpp"
#include "adjoint.hpp"
#include "logging.hpp"
#include "regularization.hpp"
#include "regularizationOptions.hpp"

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

SensitivityAndCost computeSensitivityAndCost(const FixedParameters& params, ConstMatrixRef permeabilities,
                                             ConstMatrixRef logPermeabilities, Rng& rng) {
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
        breakthroughHappened = stepForwardAndAdjointProblem(params, permeabilities, currentTimeLevel, simulationState,
                                                            randomWalks, antitheticRandomWalks, rng);
        if (enableAntitheticSampling) {
            log()->info("Antithetic sampling enabled");
        } else {
            log()->info("Antithetic sampling disabled");
        }

        if (true) {
            removeAbsorbedStates(randomWalks, numberOfRandomWalksPerParameter, sensitivities);
        }

        const Real contributionToCost = computeContributionToCost(params, simulationState);
        log()->info("contribution to cost = {}", contributionToCost);
        cost += contributionToCost;
        ++currentTimeLevel;
    } while (!breakthroughHappened && simulationState.time < params.finalTime && currentTimeLevel < params.maxNumberOfTimesteps);



    using RandomWalkIterator = typename std::list<RandomWalkState>::iterator;
    RandomWalkIterator randomWalkIterator = randomWalks.begin();
    RandomWalkIterator antitheticRandomWalkIterator = antitheticRandomWalks.begin();
    const RandomWalkIterator randomWalkEnd = randomWalks.end();


    while (randomWalkIterator != randomWalkEnd) {
        sensitivities(randomWalkIterator->parameterIndex) += randomWalkIterator->D;
        ++numberOfRandomWalksPerParameter(randomWalkIterator->parameterIndex);
        ++randomWalkIterator;
        if (enableAntitheticSampling) {
            sensitivitiesAntithetic(antitheticRandomWalkIterator->parameterIndex) += antitheticRandomWalkIterator->D;
            ++numberOfRandomWalksPerParameterAntithetic(antitheticRandomWalkIterator->parameterIndex);
            ++antitheticRandomWalkIterator;
        }
    }


    sensitivities.array() /= numberOfRandomWalksPerParameter.array().cwiseMax(1).cast<long double>();

    if (enableAntitheticSampling) {
        sensitivitiesAntithetic.array() /= numberOfRandomWalksPerParameterAntithetic.array().cwiseMax(1).cast<long double>();

    }


    Real antitheticPart = 0;
    if (enableAntitheticSampling) {
        antitheticPart = 0;
    }

    SensitivityAndCost sensitivityAndCost;
    sensitivityAndCost.sensitivity.resizeLike(sensitivities);

    sensitivityAndCost.cost = cost;

    log()->info("sensitivities norm = {}", sensitivities.norm());


    if (enableAntitheticSampling) {
        sensitivityAndCost.sensitivity =
              ((1.0 - antitheticPart) * sensitivities + antitheticPart * sensitivitiesAntithetic).cast<Real>();

        log()->info("antithetic sensitivities norm = {}", sensitivitiesAntithetic.norm());
    } else {
        sensitivityAndCost.sensitivity = sensitivities.cast<Real>();
    }

    sensitivityAndCost.sensitivity *= -1;

    if (enableRegularization) {
        log()->info("Regularization enabled");

        const Vector regularizationPenalty = deriveRegularizationPenaltyByLogPermeabilities(logPermeabilities, params.meshWidth);
        log()->info("Regularization penalty norm = {}", regularizationPenalty.norm());
        sensitivityAndCost.sensitivity += regularizationPenalty;

        const Real regularizationPenaltyCost = computeRegularizationPenalty(logPermeabilities, params.meshWidth);
        log()->info("Regularization penalty cost = {}", regularizationPenaltyCost);
        sensitivityAndCost.cost += regularizationPenaltyCost;


    } else {
        log()->info("Regularization disabled");
    }

    log()->debug("Sensitivities =\n{}", sensitivityAndCost.sensitivity);

    return sensitivityAndCost;
}


