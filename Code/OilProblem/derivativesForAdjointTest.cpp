//
// Created by Stefano Weidmann on 13.05.18.
//

#include "derivativesForAdjointFD.hpp"
#include "derivativesForAdjoint.hpp"
#include "pressure.hpp"
#include <stefCommonHeaders/logging.hpp>
#include "logging.hpp"

int main() {
    const int n = 5;

    const auto logger = stefCommonHeaders::setUpLog(spdlog::level::level_enum::debug, true, "fd.log");


    const Matrix pressures(Matrix::Random(n, n));
    pressures.setRandom();

    const Matrix saturationsWater(Matrix::Random(n, n));

    const Matrix permeabilities(Matrix::Ones(n, n));
    const Matrix logPermeabilities = permeabilities.array().log().matrix();


    const Real dynamicViscosityWater = 1;
    const Real dynamicViscosityOil = 1;


    FixedParameters params;
    params.meshWidth = 1e-4;
    params.dynamicViscosityWater = dynamicViscosityWater;
    params.dynamicViscosityOil = dynamicViscosityOil;
    params.initialPermeabilities = permeabilities;
    params.maxNumberOfTimesteps = 10;
    params.initialSaturationsWater = saturationsWater;
    params.porosity = 0.5;
    params.finalTime = 0.001;
    params.overPressureDrill = [] (const Real time) {
        return 1;
    };

    params.inflowPerUnitDepthWater = [] (const Real time) {
        return 1;
    };

    const Matrix totalMobilities = computeTotalMobilities(dynamicViscosityOil, dynamicViscosityWater);
    const Matrix totalMobilitiesBySatWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, saturationsWater, dynamicViscosityOil, dynamicViscosityWater);
    const SparseMatrix analyticPressureBySaturation = computePressureResidualsDerivedBySaturationWater(pressures, totalMobilities, totalMobilitiesBySatWater);

    const Matrix numericPressureBySaturation = computePressureResidualsDerivedFD(pressures, saturationsWater, logPermeabilities, Shift::ShiftWhere::SATURATIONS, params);


    log()->debug("numeric pressure by saturation =\n{}", numericPressureBySaturation);
    log()->debug("analytic pressure by saturation =\n{}", analyticPressureBySaturation);

    spdlog::drop_all();

    return 0;
}