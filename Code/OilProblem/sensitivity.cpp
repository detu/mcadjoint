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
        brokeThrough = stepForwardAndAdjointProblemTraditional(params, params.initialPermeabilities, currentTimelevel, simulationState, adjointMatrix, adjointRhs, completeC);

        cost += computeContributionToCost(params, simulationState);
        log()->info("time = {}", simulationState.time);
        if (!brokeThrough && (currentTimelevel % 10 == 0)) {

            dumpThis("adjointMatrixTrad", adjointMatrix);
            dumpThis("adjointRhsTrad", adjointRhs);
            dumpThis("completeCTrad", completeC);
            writeToMatFile();
            ASSERT(allFinite(adjointMatrix));
            ASSERT(allFinite(adjointRhs));
        }
    }
    log()->debug("broke through? {}", brokeThrough);

    const Vector adjoint = adjointMatrix.lu().solve(adjointRhs);
    dumpThis("adjointTrad", adjoint);

    const Vector sensitivities = -completeC.transpose() * adjoint;


    SensitivityAndCost sensitivityAndCost;
    sensitivityAndCost.sensitivity.resizeLike(sensitivities);

    sensitivityAndCost.cost = cost;
    sensitivityAndCost.sensitivity = sensitivities;
    applyRegularizationIfEnabled(params, logPermeabilities, sensitivityAndCost);

    dumpThis("sensitivityTrad", sensitivityAndCost.sensitivity);
    writeToMatFile();

    return sensitivityAndCost;


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
    sensitivities.array() *= currentTimeLevel;

    if (params.enableAntitheticSampling) {
        sensitivitiesAntithetic.array() /= numberOfRandomWalksPerParameterAntithetic.array().cwiseMax(1).cast<WiderReal>();
        sensitivitiesAntithetic.array() *= currentTimeLevel;
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

    applyRegularizationIfEnabled(params, logPermeabilities, sensitivityAndCost);

    log()->debug("Sensitivities =\n{}", sensitivityAndCost.sensitivity);

    return sensitivityAndCost;
}


