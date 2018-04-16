//
// Created by Stefano Weidmann on 02.04.18.
//

#include "oilProblem.hpp"
#include "fixedParameters.hpp"
#include "logging.hpp"

bool stepForwardProblem(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& currentState, VectorRef pressureRhs) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, currentState.saturationsWater);

    //LOGGER->debug("total mobilities {}", totalMobilities);
    const Real pressureDrillNow = params.overPressureDrill(currentState.time);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);

    //LOGGER->debug("pressure system {}", pressureSystem);


    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(currentState.time));

    adaptRhsForPressure(sourceAtDrillNow, pressureRhs, currentState.saturationsWater.rows(), currentState.saturationsWater.cols());

    //LOGGER->debug("pressure rhs {}", pressureRhs);
    currentState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);

    //LOGGER->debug("pressures {}", currentState.pressures.vec);
    return advanceSaturationsInTime(params, currentState.saturationsWater, currentState.pressures.map, totalMobilities, currentState.time);
}

Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater) {
    return permeabilities.array() * (saturationsWater.array().square() / dynamicViscosityWater + (1.0 - saturationsWater.array()).square() / dynamicViscosityOil);
}

