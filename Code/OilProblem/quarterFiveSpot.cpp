//
// Created by Stefano Weidmann on 10.04.18.
//

#include "oilproblem.hpp"
#include "logging.hpp"
#include <EigenSimplematio.hpp>
#include "specialCells.hpp"
#include <signal.h>
#include <cstdlib>

const int n = 100;
SimulationState simulationState(n, n);


void writeToFile(int ignored) {
    SMIO::EigenMatFile matFile("fieldsQuarterFiveSpot.mat");
    matFile.writeVariable("satWater", clamp(simulationState.saturationsWater, 0, 1));
    matFile.writeVariable("pressure", Matrix(simulationState.pressures.map));
    matFile.close();
    std::exit(1);
}

int main() {



    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);


    const Real meshWidth = 1.0 / n;

    const CellIndex wellCell = findWellCell(n, n);
    const CellIndex drillCell = findDrillCell(n, n);

    const Matrix permeabilities = Matrix::Ones(n, n);

    FixedParameters params;

    params.dynamicViscosityOil = 0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 0.0010518; // Water 20°C
    params.finalTime = 2;
    const Real atmosphericPressure = 1e5;
    params.pressureWell = [=] (const Real time) {
        return atmosphericPressure;
    };

    params.pressureDrill = [=] (const Real time) {
        return 2*atmosphericPressure;
    };



    params.porosity = 1;




    simulationState.saturationsWater.setZero();
    params.meshWidth = meshWidth;


    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 3;
    };


    Vector pressureRhs(n*n);
    pressureRhs.setZero();

    signal(SIGINT, writeToFile);
    signal(SIGTERM, writeToFile);


    while (simulationState.time < params.finalTime) {
        stepForwardProblem(params, permeabilities, simulationState, pressureRhs);
        LOGGER->debug("saturations water =\n{}", simulationState.saturationsWater);
        LOGGER->info("time = {}", simulationState.time);
    }

    LOGGER->debug("dim sat water = ({}, {})", simulationState.saturationsWater.rows(), simulationState.saturationsWater.cols());

    writeToFile(0);
    return 0;
}