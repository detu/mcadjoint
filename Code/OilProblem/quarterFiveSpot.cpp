//
// Created by Stefano Weidmann on 10.04.18.
//

#include "oilproblem.hpp"
#include "logging.hpp"
#include <EigenSimplematio.hpp>

int main() {

    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);


    const int n = 10;
    const Real meshWidth = 1.0 / n;

    const Matrix permeabilities = Matrix::Ones(n, n);

    FixedParameters params;

    params.dynamicViscosityOil = 0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 0.0010518; // Water 20°C
    params.finalTime = 0.01;
    const Real atmosphericPressure = 1e5;
    params.pressureWell = [=] (const Real time) {
        return atmosphericPressure;
    };

    params.pressureDrill = [=] (const Real time) {
        return 2*atmosphericPressure;
    };



    params.porosity = 1;



    SimulationState simulationState(n, n);

    simulationState.saturationsWater.setZero();

    params.meshWidth = meshWidth;
    CellIndex wellCell = {0, n-1};
    params.outflowPerUnitDepthWater = [&] (const Real time) {

        return 3*simulationState.saturationsWater();
    };

    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 3;
    };

    params.outflowPerUnitDepthOil = [&] (const Real time) {
        return params.inflowPerUnitDepthWater(time) - params.outflowPerUnitDepthWater(time);
    };

    Vector pressureRhs(n*n);
    pressureRhs.setZero();

    while (simulationState.time < params.finalTime) {
        stepForwardProblem(params, permeabilities, simulationState, pressureRhs);
        LOGGER->debug("saturations water =\n{}", simulationState.saturationsWater);
        LOGGER->info("time = {}", simulationState.time);
    }

    LOGGER->debug("dim sat water = ({}, {})", simulationState.saturationsWater.rows(), simulationState.saturationsWater.cols());

    SMIO::EigenMatFile matFile("satWater.mat");
    matFile.writeVariable("satWater", clamp(simulationState.saturationsWater, 0, 1));
    matFile.writeVariable("pressure", Matrix(simulationState.pressures.map));

    return 0;
}