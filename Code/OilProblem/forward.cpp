//
// Created by Stefano Weidmann on 02.04.18.
//

#include "oilProblem.hpp"
#include "fixedParameters.hpp"
#include "logging.hpp"

void stepForwardProblem(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& currentState, VectorRef pressureRhs) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, currentState.saturationsWater);

    LOGGER->debug("total mobilities {}", totalMobilities);
    const Real pressureDrillNow = params.pressureDrill(currentState.time);
    const Real pressureWellNow = params.pressureWell(currentState.time);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);

    LOGGER->debug("pressure system {}", pressureSystem);

    const Real waterSourceAtWellNow = -std::abs(params.outflowPerUnitDepthWater(currentState.time));
    const Real oilSourceAtWellNow = -std::abs(params.outflowPerUnitDepthOil(currentState.time));
    const Real totalSourceAtWellNow = waterSourceAtWellNow + oilSourceAtWellNow;

    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(currentState.time));
    adaptRhsForPressure(totalSourceAtWellNow, sourceAtDrillNow, pressureRhs, currentState.saturationsWater.rows(), currentState.saturationsWater.cols());

    LOGGER->debug("pressure rhs {}", pressureRhs);
    currentState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs, currentState.pressures.vec);

    LOGGER->debug("pressures {}", currentState.pressures.vec);
    advanceSaturationsInTime(params, currentState.saturationsWater, currentState.pressures.map, totalMobilities, currentState.time);
}

Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater) {
    return permeabilities.array() * (saturationsWater.array().square() / dynamicViscosityWater + (1.0 - saturationsWater.array()).square() / dynamicViscosityOil);
}

