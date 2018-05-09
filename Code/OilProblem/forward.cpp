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

#ifdef MULTITHREADED
#include <omp.h>
#endif


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


bool stepForwardAndAdjointProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                  const int currentTimelevel, SimulationState& simulationState,
                                  std::vector<RandomWalkState>& randomWalks, std::vector<Rng>& rngs) {
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

    const Real pressureDrillNow = params.overPressureDrill(simulationState.time);
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


    const Real firstTimestep = getFirstTimestep();
    const SparseMatrix pressureResidualsByLogPermeabilities = computePressureResidualsByLogPermeability(simulationState.pressures.map, totalMobilities);
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


    const SparseMatrix saturationResidualsByLogPermeabilities = computeSaturationsWaterResidualsByLogPermeability(pressureDerivativesX, pressureDerivativesY, totalMobilities, firstTimestep, params.meshWidth);
    dumpThis("saturationResidualsByLogPermeabilities", saturationResidualsByLogPermeabilities);



    SparseMatrix pressureResidualsByPressures = computePressureResidualsDerivedByPressure(pressureSystem);


    const Matrix totalMobilitiesDerivedBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures.map, totalMobilities, totalMobilitiesDerivedBySaturationsWater);

    SparseMatrix saturationsWaterResidualsByPressures = computeSaturationWaterResidualsDerivedByPressure(pressureSystem, fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, timestep, params.meshWidth);
    const Matrix fluxFunctionFactorDerivatives = computeFluxFunctionFactorDerivatives(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    SparseMatrix saturationsWaterResidualsBySaturationsWater = computeSaturationWaterResidualsDerivedBySaturationWater(
          fluxFunctionFactors, fluxFunctionFactorDerivatives,
          darcyVelocitiesX, darcyVelocitiesY,
          pressureDerivativesX, pressureDerivativesY,
          totalMobilities, totalMobilitiesDerivedBySaturationsWater,
          timestep, params.meshWidth);

    const DiagonalMatrix inverseDiagonalPressureByPressure = extractInverseDiagonalMatrix(pressureResidualsByPressures);


    const int drillCellLinearIndex = drillCell.linearIndex(numberOfRows);
    Real correspondingFactorForRhs = inverseDiagonalPressureByPressure.diagonal()(drillCellLinearIndex);






    //const DiagonalMatrix inverseDiagonalSaturationBySaturation = extractInverseDiagonalMatrix(saturationsWaterResidualsBySaturationsWater);


    SparseMatrix correctedPressureResidualsByPressures = pressureResidualsByPressures * inverseDiagonalPressureByPressure ;
    SparseMatrix correctedPressureResidualsBySaturationsWater = pressureResidualsBySaturationsWater;
    SparseMatrix correctedSaturationsWaterResidualsByPressures = saturationsWaterResidualsByPressures * inverseDiagonalPressureByPressure;
    SparseMatrix correctedSaturationsWaterResidualsBySaturationsWater = saturationsWaterResidualsBySaturationsWater;


    constexpr bool useConvergenceFactor = true;
    Real convergenceFactor = 1;
    if (useConvergenceFactor) {

        convergenceFactor = 1.0 / (sumOfAbsEntries(correctedPressureResidualsByPressures)); // + sumOfAbsEntries(correctedPressureResidualsBySaturationsWater)
                                        //+ sumOfAbsEntries(correctedSaturationsWaterResidualsByPressures) + sumOfAbsEntries(correctedSaturationsWaterResidualsBySaturationsWater) + numberOfParameters);

    }



    dumpThis("correctedPressureResidualsByPressures", correctedPressureResidualsByPressures);
    dumpThis("correctedPressureResidualsBySaturationsWater", correctedPressureResidualsBySaturationsWater);
    dumpThis("correctedSaturationsWaterResidualsByPressures", correctedSaturationsWaterResidualsByPressures);
    dumpThis("correctedSaturationsWaterResidualsBySaturationsWater",
             correctedSaturationsWaterResidualsBySaturationsWater);



    const BVectorSurrogate b(computedPressureAtDrillCell * correspondingFactorForRhs, measuredPressureAtDrillCell * correspondingFactorForRhs, numberOfRows, numberOfCols);
    const CMatrixSurrogate c(pressureResidualsByLogPermeabilities, saturationResidualsByLogPermeabilities,
                             numberOfRows, numberOfCols);

    constexpr bool startAddingRandomWalksAtBeginning = true;
    const bool shouldAddRandomWalks = startAddingRandomWalksAtBeginning || (simulationState.saturationsWater.diagonal(-1).array() > 0).any();

    if (shouldAddRandomWalks) {
        addNewRandomWalks(numberOfRows, numberOfCols, numberOfParameters, currentTimelevel, b, c, randomWalks,
                          rngs[0]);
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
    log()->debug("convergence factor = {}", convergenceFactor);

    #ifdef MULTITHREADED
    #pragma omp parallel
    #endif
    {
        #ifdef MULTITHREADED
        Rng& rng = rngs[omp_get_thread_num()];
        #else
        Rng& rng = rngs[0];
        #endif

        #ifdef MULTITHREADED
        #pragma omp for schedule(dynamic) reduction(+: advancedRandomWalks)
        #endif
        for (int randomWalkIndex = 0; randomWalkIndex < randomWalks.size(); ++randomWalkIndex) {
            RandomWalkState& randomWalk = randomWalks[randomWalkIndex];



            bool stillInTheSameTimestep = false;
            do {
                stillInTheSameTimestep = transitionState(randomWalk, b, correctedPressureResidualsByPressures,
                                                         correctedPressureResidualsBySaturationsWater,
                                                         correctedSaturationsWaterResidualsByPressures,
                                                         correctedSaturationsWaterResidualsBySaturationsWater,
                                                         convergenceFactor,
                                                         numberOfRows, numberOfCols, rng);
            } while (stillInTheSameTimestep);

            ++advancedRandomWalks;

            if (outputProgressTransitioning) {
                log()->info("Random walk {} / {}", advancedRandomWalks, randomWalks.size());
            }
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

    const Real pressureDrillNow = params.overPressureDrill(simulationState.time);
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


    const Real firstTimestep = getFirstTimestep();
    const SparseMatrix pressureResidualsByLogPermeabilities = computePressureResidualsByLogPermeability(simulationState.pressures.map, totalMobilities);
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


    const SparseMatrix saturationResidualsByLogPermeabilities = computeSaturationsWaterResidualsByLogPermeability(pressureDerivativesX, pressureDerivativesY, totalMobilities, firstTimestep, params.meshWidth);
    dumpThis("saturationResidualsByLogPermeabilities", saturationResidualsByLogPermeabilities);



    SparseMatrix pressureResidualsByPressures = computePressureResidualsDerivedByPressure(pressureSystem);


    const Matrix totalMobilitiesDerivedBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures.map, totalMobilities, totalMobilitiesDerivedBySaturationsWater);

    SparseMatrix saturationsWaterResidualsByPressures = computeSaturationWaterResidualsDerivedByPressure(pressureSystem, fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, timestep, params.meshWidth);
    const Matrix fluxFunctionFactorDerivatives = computeFluxFunctionFactorDerivatives(simulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    SparseMatrix saturationsWaterResidualsBySaturationsWater = computeSaturationWaterResidualsDerivedBySaturationWater(
          fluxFunctionFactors, fluxFunctionFactorDerivatives,
          darcyVelocitiesX, darcyVelocitiesY,
          pressureDerivativesX, pressureDerivativesY,
          totalMobilities, totalMobilitiesDerivedBySaturationsWater,
          timestep, params.meshWidth);

    const int n = numberOfRows;
    const int stateSize = 2 * n * n;
    auto adjointMatrixDiagonalBlock = adjointMatrix.block(stateSize * currentTimelevel, stateSize * currentTimelevel, stateSize, stateSize);
    adjointMatrixDiagonalBlock << Matrix(pressureResidualsByPressures.transpose()), Matrix(saturationsWaterResidualsByPressures.transpose()),
          Matrix::Zero(n*n, n*n), Matrix::Zero(n*n, n*n);



    if (stateSize * (currentTimelevel+1) < adjointMatrix.cols()) {
        auto adjointMatrixOffDiagonalBlock = adjointMatrix.block(stateSize * currentTimelevel,
                                                                 stateSize * (currentTimelevel + 1), stateSize,
                                                                 stateSize);

        adjointMatrixOffDiagonalBlock << Matrix::Zero(n * n, n * n), Matrix::Zero(n * n, n * n),
              Matrix(pressureResidualsBySaturationsWater.transpose()), Matrix(
              saturationsWaterResidualsBySaturationsWater.transpose());
    }




    const int drillCellLinearIndex = drillCell.linearIndex(numberOfRows);

    auto adjointRhsSegment = adjointRhs.segment(stateSize * currentTimelevel, stateSize);
    const BVectorSurrogate b(computedPressureAtDrillCell, measuredPressureAtDrillCell, numberOfRows, numberOfCols);

    for (int i = 0; i < stateSize; ++i) {
        adjointRhsSegment(i) = b(i);
    }



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