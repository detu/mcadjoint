#pragma once
#include "typedefs.hpp"
#include "cellindex.hpp"

SparseMatrix computePressureResidualsDerivedByPressure(const SparseMatrix& pressureSystem);
Matrix computeTotalMobilitiesDerivedBySaturationsWater(ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const Real dynamicViscosityOil, const Real dynamicViscosityWater);
SparseMatrix computePressureResidualsDerivedBySaturationWater(ConstMatrixRef pressures, ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater);
Matrix computeFluxFunctionFactorDerivatives(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);
SparseMatrix computeSaturationWaterResidualsDerivedBySaturationWater(
      ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxFunctionDerivatives,
      ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
      ConstMatrixRef pressureDerivativesX, ConstMatrixRef pressureDerivativesY,
      ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater,
      const Real timestep, const Real meshWidth
);


SparseMatrix computeSaturationsWaterResidualsDerivedByPressure(const SparseMatrix& pressureSystem,
                                                               ConstMatrixRef fluxFunctionFactors,
                                                               ConstMatrixRef darcyVelocitiesX,
                                                               ConstMatrixRef darcyVelocitiesY,
                                                               ConstMatrixRef mobilities,
                                                               const Real timestep, const Real meshWidth);
SparseMatrix computePressureResidualsDerivedByLogPermeability(ConstMatrixRef pressures, ConstMatrixRef totalMobilities);
SparseMatrix computeSaturationsWaterResidualsDerivedByLogPermeability(ConstMatrixRef pressureGradientsX,
                                                                      ConstMatrixRef pressureGradientsY,
                                                                      ConstMatrixRef darcyVelocitiesX,
                                                                      ConstMatrixRef darcyVelocitiesY,
                                                                      ConstMatrixRef mobilities,
                                                                      ConstMatrixRef fluxFunctionFactors,
                                                                      const Real timestep, const Real meshWidth);

inline static Real hmeanDerivedBySecond(const Real a, const Real b) {
    return 2*std::pow(a/(a+b), 2);
}