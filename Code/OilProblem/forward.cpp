//
// Created by Stefano Weidmann on 02.04.18.
//

#include "oilProblem.hpp"
#include "fixedParameters.hpp"

#ifdef TODO
SimulationState stepForwardProblem(const FixedParameters& params, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const Real time) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, saturationsWater);

    SimulationState newState;

    const Real pressureDrillNow = params.pressureDrill(time);
    const Real pressureWellNow = params.pressureWell(time);
    const SparseMatrix transmissibilities = assembleAugmentedTransmissibilityMatrix(totalMobilities);
    const Vector totalSources = -assemblePressureSourceVector(pressureDrillNow, pressureWellNow);

    newState.pressures = solvePressurePoissonProblem(transmissibilities, totalSources);
}
    #endif

Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater) {
    return permeabilities.array() * (saturationsWater.array().square() / dynamicViscosityWater + (1.0 - saturationsWater.array()).square() / dynamicViscosityOil);
}

