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


bool stepForwardProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                        SimulationState& currentState) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, currentState.saturationsWater);

    //log()->debug("total mobilities {}", totalMobilities);
    const Real pressureDrillNow = params.overPressureDrill(currentState.time);
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
    b(numberOfRows) = 2 * (computedPressureAtDrillCell - measuredPressureAtDrillCell);
    return b;
}

bool stepForwardAndAdjointProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                  const int currentTimelevel, SimulationState& simulationState,
                                  std::vector<RandomWalkState>& randomWalks, Rng& rng) {
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


    Eigen::SparseLU<SparseMatrix> pressureSolver(pressureSystem);

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









    //const DiagonalMatrix inverseDiagonalSaturationBySaturation = extractInverseDiagonalMatrix(saturationsWaterResidualsBySaturationsWater);


    SparseMatrix correctedPressureResidualsByPressures(numberOfParameters, numberOfParameters);
    correctedPressureResidualsByPressures.setIdentity();

    SparseMatrix correctedPressureResidualsBySaturationsWater = pressureResidualsBySaturationsWater;
    const SparseMatrix saturationsWaterResidualsByPressuresTransposed = saturationsWaterResidualsByPressures.transpose();
    SparseMatrix correctedSaturationsWaterResidualsByPressures = pressureSolver.solve(saturationsWaterResidualsByPressuresTransposed).transpose();
    SparseMatrix correctedSaturationsWaterResidualsBySaturationsWater = saturationsWaterResidualsBySaturationsWater;





    dumpThis("correctedPressureResidualsByPressures", correctedPressureResidualsByPressures);
    dumpThis("correctedPressureResidualsBySaturationsWater", correctedPressureResidualsBySaturationsWater);
    dumpThis("correctedSaturationsWaterResidualsByPressures", correctedSaturationsWaterResidualsByPressures);
    dumpThis("correctedSaturationsWaterResidualsBySaturationsWater",
             correctedSaturationsWaterResidualsBySaturationsWater);


    Vector b = computeBVector(computedPressureAtDrillCell, measuredPressureAtDrillCell, numberOfParameters,
                              numberOfRows, numberOfCols);

    b.head(numberOfParameters) = pressureSolver.solve(b.head(numberOfParameters));


    #ifdef JUST_COMPUTE_ADJOINT
    SparseMatrix c(numberOfParameters*2, numberOfParameters);
    for (int j = 0; j < std::min(c.cols(), c.rows()); ++j) {
            c.coeffRef(j, j) = -1;
    }

    #else
    const SparseMatrix c = concatVertically(pressureResidualsByLogPermeabilities, saturationResidualsByLogPermeabilities);
    #endif

    dumpThis("c", c);

    constexpr bool startAddingRandomWalksAtBeginning = true;
    const bool shouldAddRandomWalks = startAddingRandomWalksAtBeginning || (simulationState.saturationsWater.diagonal(-1).array() > 0).any();

    if (shouldAddRandomWalks) {
        addNewRandomWalks(numberOfRows, numberOfCols, numberOfParameters, currentTimelevel, b, c, randomWalks,
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


        for (int randomWalkIndex = 0; randomWalkIndex < randomWalks.size(); ++randomWalkIndex) {
            RandomWalkState& randomWalk = randomWalks[randomWalkIndex];



            bool stillInTheSameTimestep = false;
            do {
                stillInTheSameTimestep = transitionState(randomWalk, b, correctedPressureResidualsByPressures,
                                                         correctedPressureResidualsBySaturationsWater,
                                                         correctedSaturationsWaterResidualsByPressures,
                                                         correctedSaturationsWaterResidualsBySaturationsWater,
                                                         numberOfRows,
                                                         numberOfCols, rng);
            } while (stillInTheSameTimestep);

            ++advancedRandomWalks;

            if (outputProgressTransitioning) {
                log()->info("Random walk {} / {}", advancedRandomWalks, randomWalks.size());
            }
        }

    logStatisticsAboutRandomWalks(randomWalks);




    const Matrix saturationsWaterDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);
    simulationState.saturationsWater -= timestep * saturationsWaterDivergences;
    dumpThis("saturationsWaterDivergences", saturationsWaterDivergences);
    simulationState.time += timestep;
    drillCell(simulationState.saturationsWater) = 1;

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const bool breakthroughHappened = std::abs(wellCell(simulationState.saturationsWater)) > 1e-16;

    if (currentTimelevel % 10 == 0) {
        writeToMatFile();
    }

    return breakthroughHappened;
}


bool stepForwardAndAdjointProblemTraditional(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                  const int currentTimelevel, SimulationState& simulationState, Matrix& adjointMatrix, Vector& adjointRhs) {

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

    dumpThis("saturationsWater", simulationState.saturationsWater);
    log()->debug("permeabilities =\n{}", permeabilities);

    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, simulationState.saturationsWater);

    // solve pressure system

    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(simulationState.time));


    const Vector pressureRhs = computeRhsForPressureSystem(sourceAtDrillNow, simulationState.saturationsWater.rows(),
                                                           simulationState.saturationsWater.cols());
    simulationState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);


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
    writeToMatFile();

    const int n = numberOfRows;
    const int stateSize = 2 * n * n;
    const int fromDiagRow = stateSize * currentTimelevel;
    const int fromDiagCol = fromDiagRow;
    log()->debug("Diagonal block ({}-{}) x ({}-{})", fromDiagRow, fromDiagRow + stateSize, fromDiagCol, fromDiagCol + stateSize);
    adjointMatrix.block(fromDiagRow, fromDiagCol, stateSize/2, stateSize/2) = pressureResidualsByPressures.transpose();
    adjointMatrix.block(fromDiagRow, fromDiagCol + stateSize/2, stateSize/2, stateSize/2) = saturationsWaterResidualsByPressures.transpose();
    adjointMatrix.block(fromDiagRow + stateSize/2, fromDiagCol + stateSize/2, stateSize/2, stateSize/2).setIdentity();


    log()->debug("Pressure by pressure = \n{}", pressureResidualsByPressures);
    log()->debug("Saturations water residuals by pressure =\n{}", saturationsWaterResidualsByPressures);

    ASSERT(allFinite(pressureResidualsByPressures));
    ASSERT(allFinite(saturationsWaterResidualsByPressures));



    const int fromOffDiagRow = stateSize * currentTimelevel;
    const int fromOffDiagCol = fromDiagRow + stateSize;
    if (fromOffDiagCol < adjointMatrix.cols()) {


        log()->debug("Off diagonal block ({}-{}) x ({}-{})", fromOffDiagRow, fromOffDiagRow + stateSize, fromOffDiagCol, fromOffDiagCol + stateSize);


        adjointMatrix.block(fromOffDiagRow + stateSize/2, fromOffDiagCol, stateSize/2, stateSize/2) = pressureResidualsBySaturationsWater.transpose();
        adjointMatrix.block(fromOffDiagRow + stateSize/2, fromOffDiagCol + stateSize/2, stateSize/2, stateSize/2) = saturationsWaterResidualsBySaturationsWater.transpose();

    }




    const Vector b = computeBVector(computedPressureAtDrillCell, measuredPressureAtDrillCell, numberOfParameters, numberOfRows, numberOfCols);
    for (int i = 0; i < stateSize; ++i) {
        adjointRhs.segment(stateSize * currentTimelevel, stateSize)(i) = b(i);
    }



    const Matrix saturationsWaterDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);
    simulationState.saturationsWater -= timestep * saturationsWaterDivergences;
    dumpThis("saturationsWaterDivergencesTrad", saturationsWaterDivergences);
    simulationState.time += timestep;
    drillCell(simulationState.saturationsWater) = 1;

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const bool breakthroughHappened = std::abs(wellCell(simulationState.saturationsWater)) > 1e-16;

    if (currentTimelevel % 10 == 0) {
        writeToMatFile();
    }

    return breakthroughHappened;
}