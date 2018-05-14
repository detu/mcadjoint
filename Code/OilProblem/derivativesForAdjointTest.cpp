//
// Created by Stefano Weidmann on 13.05.18.
//

#include "derivativesForAdjointFD.hpp"
#include "derivativesForAdjoint.hpp"
#include "pressure.hpp"
#include "darcyVelocity.hpp"
#include "saturation.hpp"
#include <stefCommonHeaders/logging.hpp>
#include "logging.hpp"

int main() {
    const int n = 2;

    auto logger = stefCommonHeaders::setUpLog<stefCommonHeaders::NoMutex>(spdlog::level::level_enum::debug, true, "fd.log");


    const Matrix pressures(Matrix::Random(n, n));

    logger->debug("pressure =\n{}", pressures);

    ASSERT(n == 2);
    Matrix saturationsWater(n, n);
    saturationsWater << 0.7, 0.1,
                        0.01, 0.8;

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

    const Real timestep = 1e-5;

    const Matrix totalMobilities = computeTotalMobilities(dynamicViscosityOil, dynamicViscosityWater, permeabilities, saturationsWater);
    logger->debug("total Mobilities =\n{}", totalMobilities);

    const Matrix totalMobilitiesBySatWater = deriveTotalMobilitiesBySaturations(permeabilities, saturationsWater, dynamicViscosityOil, dynamicViscosityWater);
    logger->debug("tot mob sat water =\n{}", totalMobilitiesBySatWater);
    const SparseMatrix analyticPressureBySaturation = derivePressureResidualsBySaturations(pressures, totalMobilities, totalMobilitiesBySatWater);

    const Matrix numericPressureBySaturation = deriveResidualsWithFiniteDifferences(pressures, saturationsWater,
                                                                                    logPermeabilities,
                                                                                    WhichResidual::PRESSURE, Shift::ShiftWhere::SATURATIONS,
                                                                                    params, timestep);

    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const SparseMatrix analyticPressureByPressure = derivePressureResidualsByPresures(pressureSystem);
    const Matrix numericPressureByPressure = deriveResidualsWithFiniteDifferences(pressures, saturationsWater,
                                                                                  logPermeabilities, WhichResidual::PRESSURE,
                                                                                  Shift::ShiftWhere::PRESSURES, params, timestep);


    const SparseMatrix analyticPressureByLogPermeabilities = derivePressureResidualsByLogPermeabilities(pressures, totalMobilities);
    const Matrix numericPressureByLogPermeabilities = deriveResidualsWithFiniteDifferences(pressures, saturationsWater,
                                                                                           logPermeabilities,
                                                                                           WhichResidual::PRESSURE, Shift::ShiftWhere::LOG_PERMEABILITIES,
                                                                                           params, timestep);

    const Matrix numericSaturationBySaturation = deriveResidualsWithFiniteDifferences(pressures, saturationsWater, logPermeabilities, WhichResidual::SATURATION, Shift::ShiftWhere::SATURATIONS, params, timestep);

    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    const Matrix fluxFunctionDerivatives = deriveFluxFunctionFactorsBySaturations(saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);
    const Matrix pressureDerivativesX = computeXDerivative(pressures, params.meshWidth);
    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix pressureDerivativesY = computeYDerivative(pressures, params.meshWidth);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);
    const Matrix analyticSaturationBySaturation = deriveSaturationResidualsBySaturations(fluxFunctionFactors, fluxFunctionDerivatives,
                                                                                         darcyVelocitiesX, darcyVelocitiesY,
                                                                                         pressureDerivativesX, pressureDerivativesY,
                                                                                         totalMobilities, totalMobilitiesBySatWater, timestep, params.meshWidth);


    const Matrix analyticSaturationByPressure = deriveSaturationResidualsByPressures(pressureSystem, fluxFunctionFactors, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, timestep, params.meshWidth);
    const Matrix numericSaturationByPressure = deriveResidualsWithFiniteDifferences(pressures, saturationsWater, logPermeabilities, WhichResidual::SATURATION, Shift::ShiftWhere::PRESSURES, params, timestep);

    const Matrix analyticSaturationByLogPermeabilities = deriveSaturationResidualsByLogPermeabilities(pressureDerivativesX, pressureDerivativesY, darcyVelocitiesX, darcyVelocitiesY, totalMobilities, fluxFunctionFactors, timestep, params.meshWidth);
    const Matrix numericSaturationByLogPermeabilities = deriveResidualsWithFiniteDifferences(pressures, saturationsWater, logPermeabilities, WhichResidual::SATURATION, Shift::ShiftWhere::LOG_PERMEABILITIES, params, timestep);

    logger->debug("numeric pressure by pressure =\n{}", numericPressureByPressure);
    logger->debug("analytic pressure by pressure =\n{}", analyticPressureByPressure);

    logger->debug("difference pressure by pressure =\n{}", numericPressureByPressure - analyticPressureByPressure);

    logger->debug("max error pressure by pressure =\n{}", Matrix(numericPressureByPressure - analyticPressureByPressure).array().abs().maxCoeff());


    logger->debug("numeric pressure by saturation =\n{}", numericPressureBySaturation);
    logger->debug("analytic pressure by saturation =\n{}", analyticPressureBySaturation);

    logger->debug("difference pressure by saturation =\n{}", numericPressureBySaturation - analyticPressureBySaturation);

    logger->debug("max error pressure by saturation =\n{}", Matrix(numericPressureBySaturation - analyticPressureBySaturation).array().abs().maxCoeff());


    logger->debug("numeric pressure by logPermeability =\n{}", numericPressureByLogPermeabilities);
    logger->debug("analytic pressure by logPermeability =\n{}", analyticPressureByLogPermeabilities);

    logger->debug("difference pressure by logPermeability =\n{}", numericPressureByLogPermeabilities - analyticPressureByLogPermeabilities);

    logger->debug("max error pressure by logPermeability =\n{}", Matrix(numericPressureByLogPermeabilities - analyticPressureByLogPermeabilities).array().abs().maxCoeff());


    logger->debug("numeric saturation by saturation =\n{}", numericSaturationBySaturation);
    logger->debug("analytic saturation by saturation =\n{}", analyticSaturationBySaturation);

    logger->debug("difference saturation by saturation =\n{}", numericSaturationBySaturation - analyticSaturationBySaturation);

    logger->debug("max error saturation by saturation =\n{}", Matrix(numericSaturationBySaturation - analyticSaturationBySaturation).array().abs().maxCoeff());


    logger->debug("numeric saturation by pressure =\n{}", numericSaturationByPressure);
    logger->debug("analytic saturation by pressure =\n{}", analyticSaturationByPressure);

    logger->debug("difference saturation by pressure =\n{}", numericSaturationByPressure - analyticSaturationByPressure);

    logger->debug("max error saturation by pressure =\n{}", Matrix(numericSaturationByPressure - analyticSaturationByPressure).array().abs().maxCoeff());


    logger->debug("numeric saturation by logPermeability =\n{}", numericSaturationByLogPermeabilities);
    logger->debug("analytic saturation by logPermeability =\n{}", analyticSaturationByLogPermeabilities);

    logger->debug("difference saturation by logPermeability =\n{}", numericSaturationByLogPermeabilities - analyticSaturationByLogPermeabilities);

    logger->debug("max error saturation by logPermeability =\n{}", Matrix(numericSaturationByLogPermeabilities - analyticSaturationByLogPermeabilities).array().abs().maxCoeff());




    logger->debug("Darcy velocities X =\n{}", darcyVelocitiesX);
    logger->debug("Darcy velocities Y =\n{}", darcyVelocitiesY);

    spdlog::drop_all();

    return 0;
}