//
// Created by Stefano Weidmann on 10.04.18.
//

#include "oilproblem.hpp"
#include <stefCommonHeaders/eigenToCsv.hpp>
#include <stefCommonHeaders/stefFenv.h>

int main() {

    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
    const Real finalTime = 0.2;


    const int n = 3;
    const Real meshWidth = 1.0 / n;

    const Matrix permeabilities = Matrix::Ones(n, n);

    FixedParameters params;

    params.dynamicViscosityOil = 0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 0.0010518; // Water 20°C
    const Real atmosphericPressure = 1e5;
    params.pressureWell = [=] (const Real time) {
        return atmosphericPressure;
    };

    params.pressureDrill = [=] (const Real time) {
        return 20 * params.pressureWell(time);
    };

    params.meshWidth = meshWidth;
    params.outflowPerUnitDepthWater = [=] (const Real time) {
        return 1;
    };

    params.inflowPerUnitDepthWater = [=] (const Real time) {
        return 10;
    };

    params.outflowPerUnitDepthOil = [=] (const Real time) {
        return 100;
    };

    params.porosity = 0.5;



    SimulationState simulationState(n, n);

    simulationState.saturationsWater.setZero();

    Vector pressureRhs(n*n);
    pressureRhs.setZero();

    while (simulationState.time < finalTime) {
        stepForwardProblem(params, permeabilities, simulationState, pressureRhs);
    }

    stefCommonHeaders::writeToCsv("satWater.csv", simulationState.saturationsWater);

    return 0;
}