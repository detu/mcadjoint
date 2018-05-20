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
#include "specialCells.hpp"
#include "preconditioning.hpp"
#include "forwardOptions.hpp"

#include <iostream>
#include <random>

bool stepForwardProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                        SimulationState& currentState) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, currentState.saturationsWater);

    //log()->debug("total mobilities {}", totalMobilities);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);

    //log()->debug("pressure system {}", pressureSystem);


    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(currentState.time));


    //log()->debug("pressure rhs {}", pressureRhs);
    const Vector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, currentState.saturationsWater.rows(), currentState.saturationsWater.cols());
    currentState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);

    //log()->debug("pressures {}", currentState.pressures.vec);
    return advanceSaturationsInTime(params, currentState.saturationsWater, currentState.pressures.map, totalMobilities, currentState.time);
}


static Vector computeBVector(const Real computedPressureAtDrillCell, const Real measuredPressureAtDrillCell,
                            const int numberOfParameters, const int numberOfRows, const int numberOfCols) {

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    const int drillCellLinearIndex = drillCell.linearIndex(numberOfRows);

    Vector b(2 * numberOfParameters);
    b.setZero();
    b(drillCellLinearIndex) = 2 * (computedPressureAtDrillCell - measuredPressureAtDrillCell);
    return b;
}

bool stepForwardAndAdjointProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                  const int currentTimelevel, const int numberOfRandomWalksToAdd, SimulationState& simulationState,
                                  std::list<RandomWalkState>& randomWalks, std::list<RandomWalkState>& antitheticRandomWalks, Rng& rng) {
    const int numberOfRows = permeabilities.rows();
    const int numberOfCols = permeabilities.cols();
    const int numberOfParameters = permeabilities.size();

    const bool isFirstTimestep = simulationState.time <= 0;

    if (isFirstTimestep) {
        simulationState.saturationsWater.resizeLike(params.initialSaturationsWater);
        simulationState.saturationsWater = params.initialSaturationsWater;
    }

    dumpThis("saturationsWater", simulationState.saturationsWater);
    log()->debug("permeabilities =\n{}", permeabilities);

    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, simulationState.saturationsWater);

    // solve pressure system

    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(simulationState.time));


    const Vector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, simulationState.saturationsWater.rows(),
                                                                 simulationState.saturationsWater.cols());


    PressureSolver pressureSolver(pressureSystem);



    simulationState.pressures = Vector(pressureSolver.solve(pressureRhs));


    dumpThis("pressures", Matrix(simulationState.pressures.map));
    log()->debug("pressures =\n{}", simulationState.pressures.map);
    log()->debug("saturations Water =\n{}", simulationState.saturationsWater);


    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    const Real computedPressureAtDrillCell = drillCell(simulationState.pressures.map);
    const Real measuredPressureAtDrillCell = params.overPressureDrill(0);


    const SparseMatrix pressureResidualsByLogPermeabilities = derivePressureResidualsByLogPermeabilities(
          simulationState.pressures.map, totalMobilities);
    dumpThis("pressureResidualsByLogPermeabilities", pressureResidualsByLogPermeabilities);

    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    dumpThis("fluxFunctionFactors", fluxFunctionFactors);

    const Matrix pressureDerivativesX = computeXDerivative(simulationState.pressures.map, params.meshWidth);
    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix fluxesX = computeFluxesX(fluxFunctionFactors, darcyVelocitiesX);

    const Matrix pressureDerivativesY = computeYDerivative(simulationState.pressures.map, params.meshWidth);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);
    const Matrix fluxesY = computeFluxesY(fluxFunctionFactors, darcyVelocitiesY);

    const Real timestep = computeTimestep(fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, params.meshWidth, params.finalTime, simulationState.time);


    const SparseMatrix saturationResidualsByLogPermeabilities = deriveSaturationResidualsByLogPermeabilities(
          pressureDerivativesX, pressureDerivativesY, darcyVelocitiesX, darcyVelocitiesY, totalMobilities,
          fluxFunctionFactors,
          timestep, params.meshWidth);
    dumpThis("saturationResidualsByLogPermeabilities", saturationResidualsByLogPermeabilities);



    SparseMatrix pressureResidualsByPressures = derivePressureResidualsByPresures(pressureSystem);


    const Matrix totalMobilitiesDerivedBySaturationsWater = deriveTotalMobilitiesBySaturations(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    SparseMatrix pressureResidualsBySaturationsWater = derivePressureResidualsBySaturations(simulationState.pressures.map, totalMobilities, totalMobilitiesDerivedBySaturationsWater);

    SparseMatrix saturationsWaterResidualsByPressures = deriveSaturationResidualsByPressures(
          pressureSystem, fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, timestep,
          params.meshWidth);
    const Matrix fluxFunctionFactorDerivatives = deriveFluxFunctionFactorsBySaturations(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    SparseMatrix saturationsWaterResidualsBySaturationsWater = deriveSaturationResidualsBySaturations(
          fluxFunctionFactors, fluxFunctionFactorDerivatives,
          darcyVelocitiesX, darcyVelocitiesY,
          pressureDerivativesX, pressureDerivativesY,
          totalMobilities, totalMobilitiesDerivedBySaturationsWater,
          timestep, params.meshWidth);






    SparseMatrix correctedPressureResidualsByPressures = pressureResidualsByPressures;

    SparseMatrix correctedPressureResidualsBySaturationsWater = pressureResidualsBySaturationsWater;
    SparseMatrix correctedSaturationsWaterResidualsByPressures = saturationsWaterResidualsByPressures;
    SparseMatrix correctedSaturationsWaterResidualsBySaturationsWater = saturationsWaterResidualsBySaturationsWater;


    Vector b = computeBVector(computedPressureAtDrillCell, measuredPressureAtDrillCell, numberOfParameters,
                              numberOfRows, numberOfCols);
    log()->debug("measure pressure at drill = {}", measuredPressureAtDrillCell);
    log()->debug("computed pressure at drill = {}", computedPressureAtDrillCell);
    log()->debug("b = {}", b);
    dumpThis("b", b);
    ASSERT(allFinite(b));
    preconditionMatrices(numberOfRows, numberOfCols, correctedPressureResidualsByPressures, correctedSaturationsWaterResidualsByPressures,
                         correctedPressureResidualsBySaturationsWater, correctedSaturationsWaterResidualsBySaturationsWater, b,
                         pressureSolver, preconditionerToUse);

    log()->debug("correctedSaturationsWaterResidualsByPressures =\n{}", correctedSaturationsWaterResidualsByPressures);

    dumpThis("correctedPressureResidualsByPressures", correctedPressureResidualsByPressures);
    dumpThis("correctedPressureResidualsBySaturationsWater", correctedPressureResidualsBySaturationsWater);
    dumpThis("correctedSaturationsWaterResidualsByPressures", correctedSaturationsWaterResidualsByPressures);
    dumpThis("correctedSaturationsWaterResidualsBySaturationsWater",
             correctedSaturationsWaterResidualsBySaturationsWater);



    #ifdef JUST_COMPUTE_ADJOINT
    SparseMatrix c(numberOfParameters*2, numberOfParameters);
    for (int j = 0; j < std::min(c.cols(), c.rows()); ++j) {
            c.coeffRef(j + numberOfParameters, j) = -1;
    }

    #else
    const SparseMatrix c = concatVertically(pressureResidualsByLogPermeabilities, saturationResidualsByLogPermeabilities);
    #endif

    dumpThis("c", c);

    const bool shouldAddRandomWalks = startAddingRandomWalksAtBeginning || (simulationState.saturationsWater.diagonal(-1).array() > 0).any();

    if (shouldAddRandomWalks) {
        addNewRandomWalks(numberOfRows, numberOfCols, numberOfParameters, currentTimelevel, numberOfRandomWalksToAdd, params.enableAntitheticSampling, b, c, randomWalks, antitheticRandomWalks,
                          rng);
    }





    constexpr bool showResidualDerivatives = true;

    if (showResidualDerivatives) {
        log()->debug("pressure residuals by pressure =\n{}", pressureResidualsByPressures);
        log()->debug("pressure residuals by satwater =\n{}", pressureResidualsBySaturationsWater);
        log()->debug("satwater residuals by satwater =\n{}", saturationsWaterResidualsBySaturationsWater);
        log()->debug("satwater residuals by pressure =\n{}", saturationsWaterResidualsByPressures);

        log()->debug("corrected pressure residuals by pressure =\n{}", correctedPressureResidualsByPressures);
        log()->debug("corrected pressure residuals by satwater =\n{}", correctedPressureResidualsBySaturationsWater);
        log()->debug("corrected sat water residuals by sat water =\n{}", correctedSaturationsWaterResidualsBySaturationsWater);
        log()->debug("corrected saturation residuals by pressure =\n{}",
                      correctedSaturationsWaterResidualsByPressures);
    }

    int advancedRandomWalks = 0;
    constexpr bool outputProgressTransitioning = false;
    Real standardUniformNumber;

    std::uniform_real_distribution<Real> standardUniformDistribution;
    auto randomWalkIterator = randomWalks.begin();
    const auto endRandomWalk = randomWalks.end();
    auto antitheticRandomWalkIterator = antitheticRandomWalks.begin();
    const auto endAntitheticRandomWalk = antitheticRandomWalks.end();

    bool regularRandomWalksFinished = randomWalkIterator == endRandomWalk;
    bool antitheticRandomWalksFinished = implies(params.enableAntitheticSampling, antitheticRandomWalkIterator == endAntitheticRandomWalk);

    for (;;) {
        standardUniformNumber = standardUniformDistribution(rng);


        if (!regularRandomWalksFinished) {
            const bool stillInTheSameTimestep = transitionState(*randomWalkIterator, b, correctedPressureResidualsByPressures,
                                                          correctedPressureResidualsBySaturationsWater,
                                                          correctedSaturationsWaterResidualsByPressures,
                                                          correctedSaturationsWaterResidualsBySaturationsWater,
                                                          numberOfRows,
                                                          numberOfCols, standardUniformNumber);
            if (!stillInTheSameTimestep) {
                 ++randomWalkIterator;
                regularRandomWalksFinished = randomWalkIterator == endRandomWalk;
            }
        }

        if (!antitheticRandomWalksFinished) {
            const bool stillInTheSameTimestepAntithetic = transitionState(*antitheticRandomWalkIterator, b,
                                                               correctedPressureResidualsByPressures,
                                                               correctedPressureResidualsBySaturationsWater,
                                                               correctedSaturationsWaterResidualsByPressures,
                                                               correctedSaturationsWaterResidualsBySaturationsWater,
                                                               numberOfRows,
                                                               numberOfCols, 1.0 - standardUniformNumber);
            if (!stillInTheSameTimestepAntithetic) {
                ++antitheticRandomWalkIterator;
                antitheticRandomWalksFinished = antitheticRandomWalkIterator == endAntitheticRandomWalk;
            }
        }

        if (antitheticRandomWalksFinished && regularRandomWalksFinished) {
            ++advancedRandomWalks;

            if (outputProgressTransitioning) {
                log()->info("Random walk {} / {}", advancedRandomWalks, randomWalks.size());
            }

            break;
        }
    }


    logStatisticsAboutRandomWalks(randomWalks);

    if (params.enableAntitheticSampling) {
        logStatisticsAboutRandomWalks(antitheticRandomWalks);
    }




    const Matrix saturationsWaterDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);
    simulationState.saturationsWater -= timestep * saturationsWaterDivergences;
    dumpThis("saturationsWaterDivergences", saturationsWaterDivergences);
    simulationState.time += timestep;
    drillCell(simulationState.saturationsWater) = 1;

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const bool breakthroughHappened = std::abs(wellCell(simulationState.saturationsWater)) > 0.9;

    if (currentTimelevel % 10 == 0) {
        writeToMatFile();
    }

    return breakthroughHappened;
}


bool stepForwardAndAdjointProblemTraditional(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                  const int currentTimelevel, SimulationState& simulationState, Matrix& adjointMatrix, Vector& adjointRhs, Matrix& completeC) {

    ASSERT(allFinite(adjointMatrix));
    ASSERT(allFinite(adjointRhs));

    const int numberOfRows = permeabilities.rows();
    const int numberOfCols = permeabilities.cols();
    const int numberOfParameters = permeabilities.size();

    const bool isFirstTimestep = simulationState.time <= 0;

    if (isFirstTimestep) {
        simulationState.saturationsWater.resizeLike(params.initialSaturationsWater);
        simulationState.saturationsWater = params.initialSaturationsWater;
    }

    dumpThis("saturationsWaterTrad", simulationState.saturationsWater);
    log()->debug("permeabilities =\n{}", permeabilities);

    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, simulationState.saturationsWater);

    // solve pressure system

    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(simulationState.time));


    const Vector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, simulationState.saturationsWater.rows(),
                                                           simulationState.saturationsWater.cols());
    simulationState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);


    dumpThis("pressuresTrad", Matrix(simulationState.pressures.map));
    log()->debug("pressures =\n{}", simulationState.pressures.map);
    log()->debug("saturations Water =\n{}", simulationState.saturationsWater);


    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    const Real computedPressureAtDrillCell = drillCell(simulationState.pressures.map);
    const Real measuredPressureAtDrillCell = params.overPressureDrill(0);




    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    dumpThis("fluxFunctionFactorsTrad", fluxFunctionFactors);

    const Matrix pressureDerivativesX = computeXDerivative(simulationState.pressures.map, params.meshWidth);
    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix fluxesX = computeFluxesX(fluxFunctionFactors, darcyVelocitiesX);

    const Matrix pressureDerivativesY = computeYDerivative(simulationState.pressures.map, params.meshWidth);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);
    const Matrix fluxesY = computeFluxesY(fluxFunctionFactors, darcyVelocitiesY);

    const Real timestep = computeTimestep(fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, params.meshWidth, params.finalTime, simulationState.time);






    SparseMatrix pressureResidualsByPressures = derivePressureResidualsByPresures(pressureSystem);

    const Matrix totalMobilitiesDerivedBySaturationsWater = deriveTotalMobilitiesBySaturations(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    SparseMatrix pressureResidualsBySaturationsWater = derivePressureResidualsBySaturations(simulationState.pressures.map, totalMobilities, totalMobilitiesDerivedBySaturationsWater);
    ASSERT(allFinite(pressureResidualsBySaturationsWater));

    SparseMatrix saturationsWaterResidualsByPressures = deriveSaturationResidualsByPressures(
          pressureSystem, fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, timestep,
          params.meshWidth);
    ASSERT(allFinite(saturationsWaterResidualsByPressures));

    const Matrix fluxFunctionFactorDerivatives = deriveFluxFunctionFactorsBySaturations(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    SparseMatrix saturationsWaterResidualsBySaturationsWater = deriveSaturationResidualsBySaturations(
          fluxFunctionFactors, fluxFunctionFactorDerivatives,
          darcyVelocitiesX, darcyVelocitiesY,
          pressureDerivativesX, pressureDerivativesY,
          totalMobilities, totalMobilitiesDerivedBySaturationsWater,
          timestep, params.meshWidth);



    dumpThis("saturationsWaterResidualsBySaturationsWaterTrad", saturationsWaterResidualsBySaturationsWater);
    dumpThis("saturationsWaterResidualsByPressuresTrad", saturationsWaterResidualsByPressures);
    dumpThis("pressureResidualsBySaturationsWaterTrad", pressureResidualsBySaturationsWater);
    dumpThis("pressureResidualsByPressuresTrad", pressureResidualsByPressures);
    //writeToMatFile();

    const int stateSize = 2 * numberOfRows * numberOfCols;
    const int fromDiagRow = stateSize * currentTimelevel;
    const int fromDiagCol = fromDiagRow;
    log()->debug("Diagonal block ({}-{}) x ({}-{})", fromDiagRow, fromDiagRow + stateSize, fromDiagCol, fromDiagCol + stateSize);
    log()->debug("Pressure by pressure = \n{}", pressureResidualsByPressures);
    log()->debug("Saturations water residuals by pressure =\n{}", saturationsWaterResidualsByPressures);
    adjointMatrix.block(fromDiagRow, fromDiagCol, stateSize/2, stateSize/2) = pressureResidualsByPressures.transpose();

    Matrix densePressureResidualsByPressure = Matrix::Zero(stateSize/2, stateSize/2);
    densePressureResidualsByPressure = pressureResidualsByPressures;

    dumpThis("densePressureResidualsByPressure", densePressureResidualsByPressure);
    //writeToMatFile();
    ASSERT(allFinite(densePressureResidualsByPressure));

    ASSERT(allFinite(adjointMatrix.block(fromDiagRow, fromDiagCol, stateSize/2, stateSize/2)));


    adjointMatrix.block(fromDiagRow, fromDiagCol + stateSize/2, stateSize/2, stateSize/2) = saturationsWaterResidualsByPressures.transpose();

    ASSERT(allFinite(adjointMatrix.block(fromDiagRow, fromDiagCol + stateSize/2, stateSize/2, stateSize/2)));

    adjointMatrix.block(fromDiagRow + stateSize/2, fromDiagCol + stateSize/2, stateSize/2, stateSize/2).setIdentity();




    ASSERT(allFinite(pressureResidualsByPressures));
    ASSERT(allFinite(saturationsWaterResidualsByPressures));



    const int fromOffDiagRow = stateSize * currentTimelevel;
    const int fromOffDiagCol = fromDiagRow + stateSize;
    if (fromOffDiagCol < adjointMatrix.cols()) {


        log()->debug("Off diagonal block ({}-{}) x ({}-{})", fromOffDiagRow, fromOffDiagRow + stateSize, fromOffDiagCol, fromOffDiagCol + stateSize);


        adjointMatrix.block(fromOffDiagRow + stateSize/2, fromOffDiagCol, stateSize/2, stateSize/2) = pressureResidualsBySaturationsWater.transpose();
        ASSERT(allFinite(adjointMatrix.block(fromOffDiagRow + stateSize/2, fromOffDiagCol, stateSize/2, stateSize/2)));
        adjointMatrix.block(fromOffDiagRow + stateSize/2, fromOffDiagCol + stateSize/2, stateSize/2, stateSize/2) = saturationsWaterResidualsBySaturationsWater.transpose();
        ASSERT(allFinite(adjointMatrix.block(fromOffDiagRow + stateSize/2, fromOffDiagCol + stateSize/2, stateSize/2, stateSize/2)));
    }




    Vector b = computeBVector(computedPressureAtDrillCell, measuredPressureAtDrillCell, numberOfParameters, numberOfRows, numberOfCols);
    dumpThis("bTrad", b);
    adjointRhs.segment(stateSize * currentTimelevel, stateSize) = b;

    const SparseMatrix pressureResidualsByLogPermeabilities = derivePressureResidualsByLogPermeabilities(
          simulationState.pressures.map, totalMobilities);
    dumpThis("pressureResidualsByLogPermeabilitiesTrad", pressureResidualsByLogPermeabilities);

    const SparseMatrix saturationResidualsByLogPermeabilities = deriveSaturationResidualsByLogPermeabilities(
          pressureDerivativesX, pressureDerivativesY, darcyVelocitiesX, darcyVelocitiesY, totalMobilities,
          fluxFunctionFactors,
          timestep, params.meshWidth);
    dumpThis("saturationResidualsByLogPermeabilitiesTrad", saturationResidualsByLogPermeabilities);


    const int fromCRow = fromOffDiagRow;
    completeC.middleRows(fromCRow, stateSize/2) = pressureResidualsByLogPermeabilities;
    completeC.middleRows(fromCRow + stateSize/2, stateSize/2) = saturationResidualsByLogPermeabilities;


    const Matrix saturationsWaterDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);
    simulationState.saturationsWater -= timestep * saturationsWaterDivergences;
    dumpThis("saturationsWaterDivergencesTrad", saturationsWaterDivergences);
    simulationState.time += timestep;
    drillCell(simulationState.saturationsWater) = 1;

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const bool breakthroughHappened = std::abs(wellCell(simulationState.saturationsWater)) > 0.9;

    if (currentTimelevel % 10 == 0) {
        writeToMatFile();
    }

    return breakthroughHappened;
}