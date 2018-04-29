//
// Created by stefanow on 4/27/18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/dev.hpp>




Matrix computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities) {
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();
    const int numberOfCells = numberOfParameters;

    // TODO
    Vector pressureRhs(Vector::Zero(numberOfCells));
    SimulationState simulationState(numberOfRows, numberOfCols);

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);




    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, initialSimulationState.saturationsWater);

    Vector pressureRhs(numberOfCells);
    pressureRhs.setZero();
    std::vector<RandomWalkState> randomWalks;
    std::vector<RandomWalkState> stuckRandomWalks;
    Rng rng;


    do {


        // solve pressure system
        const Real pressureDrillNow = params.overPressureDrill(simulationState.time);
        const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
        const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(simulationState.time));
        adaptRhsForPressure(sourceAtDrillNow, pressureRhs, simulationState.saturationsWater.rows(), simulationState.saturationsWater.cols());
        simulationState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);

        const Real computedPressureAtDrillCell = drillCell(simulationState.pressures);
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

        const SparseMatrix saturationResidualsByLogPermeabilities = computeSaturationsWaterResidualsByLogPermeability(fluxesX, fluxesY, totalMobilities, firstTimestep, params.meshWidth);

        if (randomWalks.empty()) {
            const CMatrixSurrogate c(pressureResidualsByLogPermeabilities, saturationResidualsByLogPermeabilities,
                                     numberOfRows, numberOfCols);
            // initialization
            randomWalks = initializeRandomWalks(numberOfRows, numberOfCols, numberOfParameters, b, c);
        }



        const SparseMatrix pressureResidualsByPressures = computePressureResidualsDerivedByPressure(pressureSystem);

        const Matrix totalMobilitiesDerivedBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
        const SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures, totalMobilities, totalMobilitiesDerivedBySaturationsWater);

        const SparseMatrix saturationsWaterResidualsByPressure = computeSaturationWaterResidualsDerivedByPressure(pressureSystem);

        const SparseMatrix saturationsWaterResidualsBySaturationsWater = computeSaturationWaterResidualsDerivedBySaturationWater(fluxFunctionFactorDerivatives, darcyVelocitiesX, darcyVelocitiesY, timestep, params.meshWidth);

        for (RandomWalkState& randomWalk: randomWalks) {
            transitionState(randomWalk, b,
                            pressureResidualsByPressures, pressureResidualsBySaturationsWater,
                            saturationsWaterResidualsByPressure, saturationsWaterResidualsBySaturationsWater,
                            numberOfRows, numberOfCols, rng);
        }
        do {

        } while (!stuckRandomWalks.empty());
    } while (!breakthroughHappened);


    const Matrix totalMobilitiesBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    const SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures, totalMobilities, totalMobilitiesBySaturationsWater);




    return sensitivity;
}



CostFunction getCostFunctionForMinimizer(const FixedParameters& params) {
    return [&] (ConstMatrixRef logPermeabilities) {
        return computeCost(params, logPermeabilities.array().exp().matrix());
    };
}