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
#include "dumpToMatFile.hpp"
#include "utils.hpp"

bool stepForwardProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                        SimulationState& currentState) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, currentState.saturationsWater);

    //LOGGER->debug("total mobilities {}", totalMobilities);
    const Real pressureDrillNow = params.overPressureDrill(currentState.time);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);

    //LOGGER->debug("pressure system {}", pressureSystem);


    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(currentState.time));


    //LOGGER->debug("pressure rhs {}", pressureRhs);
    const Vector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, currentState.saturationsWater.rows(), currentState.saturationsWater.cols());
    currentState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);

    //LOGGER->debug("pressures {}", currentState.pressures.vec);
    return advanceSaturationsInTime(params, currentState.saturationsWater, currentState.pressures.map, totalMobilities, currentState.time);
}



bool stepForwardAndAdjointProblem(const FixedParameters& params, ConstMatrixRef permeabilities, const int currentTimelevel,
                                  SimulationState& simulationState, std::vector<RandomWalkState>& randomWalks, Rng& rng) {
    const int numberOfRows = permeabilities.rows();
    const int numberOfCols = permeabilities.cols();
    const int numberOfParameters = permeabilities.size();

    const bool isFirstTimestep = simulationState.time <= 0 || randomWalks.empty();

    if (isFirstTimestep) {
        simulationState.saturationsWater.resizeLike(params.initialSaturationsWater);
        simulationState.saturationsWater = params.initialSaturationsWater;
    }

    dumpThisOnExit("saturationsWater", simulationState.saturationsWater);
    LOGGER->debug("permeabilities =\n{}", permeabilities);

    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, simulationState.saturationsWater);

    // solve pressure system

    const Real pressureDrillNow = params.overPressureDrill(simulationState.time);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(simulationState.time));


    const Vector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, simulationState.saturationsWater.rows(),
                                                                 simulationState.saturationsWater.cols());
    simulationState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);


    dumpThisOnExit("pressures", Matrix(simulationState.pressures.map));
    LOGGER->debug("pressures =\n{}", simulationState.pressures.map);
    LOGGER->debug("saturations Water =\n{}", simulationState.saturationsWater);


    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    const Real computedPressureAtDrillCell = drillCell(simulationState.pressures.map);
    const Real measuredPressureAtDrillCell = params.overPressureDrill(0);


    const Real firstTimestep = getFirstTimestep();
    const SparseMatrix pressureResidualsByLogPermeabilities = computePressureResidualsByLogPermeability(simulationState.pressures.map, totalMobilities);

    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);

    const Matrix pressureDerivativesX = computeXDerivative(simulationState.pressures.map, params.meshWidth);
    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix fluxesX = computeFluxesX(fluxFunctionFactors, darcyVelocitiesX);

    const Matrix pressureDerivativesY = computeYDerivative(simulationState.pressures.map, params.meshWidth);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);
    const Matrix fluxesY = computeFluxesY(fluxFunctionFactors, darcyVelocitiesY);

    const Real timestep = computeTimestep(fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, params.meshWidth, params.finalTime, simulationState.time);


    const SparseMatrix saturationResidualsByLogPermeabilities = computeSaturationsWaterResidualsByLogPermeability(fluxesX, fluxesY, totalMobilities, firstTimestep, params.meshWidth);




    const SparseMatrix pressureResidualsByPressures = computePressureResidualsDerivedByPressure(pressureSystem);


    const Matrix totalMobilitiesDerivedBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    const SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures.map, totalMobilities, totalMobilitiesDerivedBySaturationsWater);

    const SparseMatrix saturationsWaterResidualsByPressures = computeSaturationWaterResidualsDerivedByPressure(pressureSystem, fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, timestep, params.meshWidth);
    const Matrix fluxFunctionFactorDerivatives = computeFluxFunctionFactorDerivatives(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    const SparseMatrix saturationsWaterResidualsBySaturationsWater = computeSaturationWaterResidualsDerivedBySaturationWater(fluxFunctionFactorDerivatives, darcyVelocitiesX, darcyVelocitiesY, timestep, params.meshWidth);


    const DiagonalMatrix inverseDiagonalPressureByPressure = extractInverseDiagonalMatrix(pressureResidualsByPressures);
    const DiagonalMatrix inverseDiagonalSaturationBySaturation = extractInverseDiagonalMatrix(saturationsWaterResidualsBySaturationsWater);

    const SparseMatrix correctedPressureResidualsByPressures = pressureResidualsByPressures * inverseDiagonalPressureByPressure;
    const SparseMatrix correctedPressureResidualsBySaturationsWater = pressureResidualsBySaturationsWater * inverseDiagonalSaturationBySaturation;
    const SparseMatrix correctedSaturationsWaterResidualsByPressures = saturationsWaterResidualsByPressures * inverseDiagonalPressureByPressure;
    const SparseMatrix correctedsaturationsWaterResidualsBySaturationsWater = saturationsWaterResidualsBySaturationsWater * inverseDiagonalSaturationBySaturation;

    dumpThisOnExit("correctedPressureResidualsByPressures", correctedPressureResidualsByPressures);
    dumpThisOnExit("correctedPressureResidualsBySaturationsWater", correctedPressureResidualsBySaturationsWater);
    dumpThisOnExit("correctedSaturationsWaterResidualsByPressures", correctedSaturationsWaterResidualsByPressures);
    dumpThisOnExit("correctedsaturationsWaterResidualsBySaturationsWater", correctedsaturationsWaterResidualsBySaturationsWater);




    const int drillCellLinearIndex = drillCell.linearIndex(numberOfRows);

    const Real targetFrobeniusNorm = 1e-6;
    const Real actualFrobeniusNorm = std::sqrt(
          frobeniusNormSquared(correctedPressureResidualsByPressures) + frobeniusNormSquared(correctedPressureResidualsBySaturationsWater)
                + frobeniusNormSquared(correctedSaturationsWaterResidualsByPressures) + frobeniusNormSquared(correctedsaturationsWaterResidualsBySaturationsWater)
    );

    const Real correspondingFactorForRhs = inverseDiagonalPressureByPressure.diagonal()(drillCellLinearIndex) * targetFrobeniusNorm / actualFrobeniusNorm;



    const BVectorSurrogate b(computedPressureAtDrillCell * correspondingFactorForRhs, measuredPressureAtDrillCell * correspondingFactorForRhs, numberOfRows, numberOfCols);
    const CMatrixSurrogate c(pressureResidualsByLogPermeabilities, saturationResidualsByLogPermeabilities,
                             numberOfRows, numberOfCols);
        // initialization
    addNewRandomWalks(numberOfRows, numberOfCols, numberOfParameters, currentTimelevel, b, c, randomWalks,
                      rng);


    constexpr bool showResidualDerivatives = true;





    if (showResidualDerivatives) {
        LOGGER->debug("pressure residuals by pressure =\n{}", pressureResidualsByPressures);
        LOGGER->debug("pressure residuals by satwater =\n{}", pressureResidualsBySaturationsWater);
        LOGGER->debug("satwater residuals by satwater =\n{}", saturationsWaterResidualsBySaturationsWater);
        LOGGER->debug("satwater residuals by pressure =\n{}", saturationsWaterResidualsByPressures);

        LOGGER->debug("corrected pressure residuals by pressure =\n{}", pressureResidualsByPressures * inverseDiagonalPressureByPressure);
        LOGGER->debug("corrected pressure residuals by satwater =\n{}", pressureResidualsBySaturationsWater * inverseDiagonalSaturationBySaturation);
        LOGGER->debug("corrected sat water residuals by sat water =\n{}", saturationsWaterResidualsBySaturationsWater * inverseDiagonalSaturationBySaturation);
        LOGGER->debug("corrected saturation residuals by pressure =\n{}",
                      saturationsWaterResidualsByPressures * inverseDiagonalPressureByPressure);
    }

    int advancedRandomWalks = 0;
    constexpr bool outputProgressTransitioning = true;
    for (RandomWalkState& randomWalk: randomWalks) {
        bool stillInTheSameTimestep = false;
        do {
            stillInTheSameTimestep = transitionState(randomWalk, b,
                                                     correctedPressureResidualsByPressures,
                                                     correctedPressureResidualsBySaturationsWater,
                                                     correctedSaturationsWaterResidualsByPressures,
                                                     correctedsaturationsWaterResidualsBySaturationsWater,
                                                     numberOfRows, numberOfCols, rng);
        } while (stillInTheSameTimestep);

        ++advancedRandomWalks;

        if (outputProgressTransitioning) {
            LOGGER->info("Random walk {} / {}", advancedRandomWalks, randomWalks.size());
        }
    }

    logStatisticsAboutRandomWalks(randomWalks);


    const Matrix saturationsWaterDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);
    simulationState.saturationsWater -= timestep * saturationsWaterDivergences;
    simulationState.time += timestep;
    drillCell(simulationState.saturationsWater) = 1;

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const bool breakthroughHappened = std::abs(wellCell(simulationState.saturationsWater)) > 1e-16;

    return breakthroughHappened;
}