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
#include "logging.hpp"
#include <cmath>

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

SensitivityAndCost computeSensitivityAndCost(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                             Rng& rng) {
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();

    SimulationState simulationState(numberOfRows, numberOfCols);
    SensitivityAndCost sensitivityAndCost = {Vector::Zero(numberOfParameters), 0};
    std::vector<RandomWalkState> randomWalks;

    bool breakthroughHappened = false;
    int currentTimeLevel = 0;
    do {
        log()->info("-----------------------------------");
        log()->info("time = {}", simulationState.time);
        breakthroughHappened = stepForwardAndAdjointProblem(params, permeabilities, currentTimeLevel, simulationState,
                                                            randomWalks, rng);
        const Real contributionToCost = computeContributionToCost(params, simulationState);
        log()->info("contribution to cost = {}", contributionToCost);
        sensitivityAndCost.cost += contributionToCost;
        ++currentTimeLevel;
    } while (!breakthroughHappened && simulationState.time < params.finalTime && currentTimeLevel < params.maxNumberOfTimesteps);

    Vector numberOfRandomWalksPerParameter(Vector::Zero(numberOfParameters));


    for (const RandomWalkState& randomWalk: randomWalks) {
        ASSERT(!std::isnan(randomWalk.D));
        sensitivityAndCost.sensitivity(randomWalk.parameterIndex) += randomWalk.D;
        ++numberOfRandomWalksPerParameter(randomWalk.parameterIndex);
    }

    log()->debug("Sensitivities before division =\n{}", sensitivityAndCost.sensitivity);

    sensitivityAndCost.sensitivity.array() /= numberOfRandomWalksPerParameter.array().cwiseMax(1);

    sensitivityAndCost.sensitivity.array() *= -1; // see equation (3)

    return sensitivityAndCost;
}


