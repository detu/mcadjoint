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
    const Real measuredPressureAtDrill = parameters.overPressureDrill(currentSimulationState.time);
    return std::pow(measuredPressureAtDrill - computedPressureAtDrill, 2);
}

SensitivityAndCost computeSensitivityAndCost(const FixedParameters& params, ConstMatrixRef permeabilities) {
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();
    const int numberOfCells = numberOfParameters;

    SimulationState simulationState(numberOfRows, numberOfCols);
    SensitivityAndCost sensitivityAndCost = {Vector::Zero(numberOfParameters), 0};
    std::vector<RandomWalkState> randomWalks;
    Rng rng;
    bool breakthroughHappened = false;
    do {
        LOGGER->info("time = {}", simulationState.time);
        breakthroughHappened = stepForwardAndAdjointProblem(params, permeabilities, simulationState, randomWalks, rng);
        sensitivityAndCost.cost += computeContributionToCost(params, simulationState);
    } while (!breakthroughHappened && simulationState.time < params.finalTime);

    Vector numberOfRandomWalksPerParameter(Vector::Zero(numberOfParameters));


    for (const RandomWalkState& randomWalk: randomWalks) {
        sensitivityAndCost.sensitivity(randomWalk.parameterIndex) += randomWalk.D;
        ++numberOfRandomWalksPerParameter(randomWalk.parameterIndex);
    }

    sensitivityAndCost.sensitivity.array() /= numberOfRandomWalksPerParameter.array();


    return sensitivityAndCost;
}


