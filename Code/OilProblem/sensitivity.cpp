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
    SensitivityAndCost sensitivityAndCost = {Vector::Zero(numberOfParameters), 0};
    std::list<RandomWalkState> randomWalks;
    std::list<RandomWalkState> antitheticRandomWalks;

    bool breakthroughHappened = false;
    int currentTimeLevel = 0;

    Eigen::VectorXi numberOfRemovedAbsorbedStates = Eigen::VectorXi::Zero(numberOfParameters);
    Vector sumOfDValuesOfAbsorbedStates = Vector::Zero(numberOfParameters);

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

        removeAbsorbedStates(randomWalks, numberOfRemovedAbsorbedStates, sumOfDValuesOfAbsorbedStates);

        const Real contributionToCost = computeContributionToCost(params, simulationState);
        log()->info("contribution to cost = {}", contributionToCost);
        sensitivityAndCost.cost += contributionToCost;
        ++currentTimeLevel;
    } while (!breakthroughHappened && simulationState.time < params.finalTime && currentTimeLevel < params.maxNumberOfTimesteps);

    Eigen::VectorXi numberOfRandomWalksPerParameter = Eigen::VectorXi(numberOfParameters);


    using RandomWalkIterator = typename std::list<RandomWalkState>::iterator;
    RandomWalkIterator randomWalkIterator = randomWalks.begin();
    RandomWalkIterator antitheticRandomWalkIterator = antitheticRandomWalks.begin();
    const RandomWalkIterator randomWalkEnd = randomWalks.end();

    const auto processRandomWalkIterator = [&sensitivityAndCost, &numberOfRandomWalksPerParameter](const RandomWalkIterator& randomWalkIterator) -> void {
          ASSERT(!std::isnan(randomWalkIterator->D));
          sensitivityAndCost.sensitivity(randomWalkIterator->parameterIndex) += randomWalkIterator->D;
          ++numberOfRandomWalksPerParameter(randomWalkIterator->parameterIndex);
    };

    while (randomWalkIterator != randomWalkEnd) {
        processRandomWalkIterator(randomWalkIterator);
        ++randomWalkIterator;
        if (enableAntitheticSampling) {
            processRandomWalkIterator(antitheticRandomWalkIterator);
            ++antitheticRandomWalkIterator;
        }
    }

    sensitivityAndCost.sensitivity += sumOfDValuesOfAbsorbedStates;
    numberOfRandomWalksPerParameter += numberOfRemovedAbsorbedStates;


    sensitivityAndCost.sensitivity.array() /= numberOfRandomWalksPerParameter.array().cwiseMax(1).cast<Real>();

    sensitivityAndCost.sensitivity.array() *= -1; // see equation (3)


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


