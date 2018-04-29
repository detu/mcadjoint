//
// Created by Stefano Weidmann on 02.04.18.
//

#include "fixedParameters.hpp"
#include "logging.hpp"
#include "pressure.hpp"
#include "forward.hpp"
#include "saturation.hpp"
#include "darcyVelocity.hpp"
#include "derivativesForAdjoint.hpp"
#include "adjoint.hpp"

bool stepForwardProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                        SimulationState& currentState) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, currentState.saturationsWater);

    //LOGGER->debug("total mobilities {}", totalMobilities);
    const Real pressureDrillNow = params.overPressureDrill(currentState.time);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);

    //LOGGER->debug("pressure system {}", pressureSystem);


    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(currentState.time));

    computeRhsForPressureSystem(sourceAtDrillNow, currentState.saturationsWater.rows(),
                                currentState.saturationsWater.cols());

    //LOGGER->debug("pressure rhs {}", pressureRhs);
    const SparseVector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, currentState.saturationsWater.rows(), currentState.saturationsWater.cols());
    currentState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);

    //LOGGER->debug("pressures {}", currentState.pressures.vec);
    return advanceSaturationsInTime(params, currentState.saturationsWater, currentState.pressures.map, totalMobilities, currentState.time);
}



bool stepForwardAndAdjointProblem(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& simulationState, std::vector<RandomWalkState>& randomWalks, Rng& rng) {
    const int numberOfRows = permeabilities.rows();
    const int numberOfCols = permeabilities.cols();
    const int numberOfParameters = permeabilities.size();

    const bool isFirstTimestep = simulationState.time <= 0 || randomWalks.empty();

    if (isFirstTimestep) {
        initializeSaturationsWater(simulationState.saturationsWater, numberOfRows, numberOfCols);
    }

    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, simulationState.saturationsWater);

    // solve pressure system

    const Real pressureDrillNow = params.overPressureDrill(simulationState.time);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(simulationState.time));


    const SparseVector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, simulationState.saturationsWater.rows(),
                                                                 simulationState.saturationsWater.cols());
    simulationState.pressures.vec = solvePressurePoissonProblem(pressureSystem, pressureRhs);

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    const Real computedPressureAtDrillCell = drillCell(simulationState.pressures.map);
    const Real measuredPressureAtDrillCell = params.overPressureDrill(0);


    const BVectorSurrogate b(computedPressureAtDrillCell, measuredPressureAtDrillCell, numberOfRows, numberOfCols);
    const Real firstTimestep = getFirstTimestep();
    const SparseMatrix pressureResidualsByLogPermeabilities = computePressureResidualsByLogPermeability(simulationState.pressures, totalMobilities);

    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);

    const Matrix pressureDerivativesX = computeXDerivative(simulationState.pressures, params.meshWidth);
    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix fluxesX = computeFluxesX(fluxFunctionFactors, darcyVelocitiesX);

    const Matrix pressureDerivativesY = computeYDerivative(simulationState.pressures, params.meshWidth);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);
    const Matrix fluxesY = computeFluxesY(fluxFunctionFactors, darcyVelocitiesY);

    const Real timestep = computeTimestep(fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, params.meshWidth, params.finalTime, simulationState.time);


    const SparseMatrix saturationResidualsByLogPermeabilities = computeSaturationsWaterResidualsByLogPermeability(fluxesX, fluxesY, totalMobilities, firstTimestep, params.meshWidth);

    if (isFirstTimestep) {
        const CMatrixSurrogate c(pressureResidualsByLogPermeabilities, saturationResidualsByLogPermeabilities,
                                 numberOfRows, numberOfCols);
        // initialization
        randomWalks = initializeRandomWalks(numberOfRows, numberOfCols, numberOfParameters, b, c);
    }



    const SparseMatrix pressureResidualsByPressures = computePressureResidualsDerivedByPressure(pressureSystem);

    const Matrix totalMobilitiesDerivedBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    const SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures, totalMobilities, totalMobilitiesDerivedBySaturationsWater);

    const SparseMatrix saturationsWaterResidualsByPressure = computeSaturationWaterResidualsDerivedByPressure(pressureSystem, fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, timestep, params.meshWidth);
    const Matrix fluxFunctionFactorDerivatives = computeFluxFunctionFactorDerivatives(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    const SparseMatrix saturationsWaterResidualsBySaturationsWater = computeSaturationWaterResidualsDerivedBySaturationWater(fluxFunctionFactorDerivatives, darcyVelocitiesX, darcyVelocitiesY, timestep, params.meshWidth);

    for (RandomWalkState& randomWalk: randomWalks) {
        bool stillInTheSameTimestep = false;
        do {
            stillInTheSameTimestep = transitionState(randomWalk, b,
                                                     pressureResidualsByPressures,
                                                     pressureResidualsBySaturationsWater,
                                                     saturationsWaterResidualsByPressure,
                                                     saturationsWaterResidualsBySaturationsWater,
                                                     numberOfRows, numberOfCols, rng);
        } while (stillInTheSameTimestep);
    }

    const Matrix saturationsWaterDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);
    simulationState.saturationsWater -= timestep * saturationsWaterDivergences;
    simulationState.time += timestep;
    drillCell(simulationState.saturationsWater) = 1;

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const bool breakthroughHappened = std::abs(wellCell(simulationState.saturationsWater)) > 1e-16;

    return breakthroughHappened;
}