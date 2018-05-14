#pragma once
#include "typedefs.hpp"
#include "cellindex.hpp"

SparseMatrix derivePressureResidualsByPresures(const SparseMatrix& pressureSystem);
Matrix deriveTotalMobilitiesBySaturations(ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const Real dynamicViscosityOil, const Real dynamicViscosityWater);
SparseMatrix derivePressureResidualsBySaturations(ConstMatrixRef pressures, ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater);
Matrix deriveFluxFunctionFactorsBySaturations(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);
SparseMatrix deriveSaturationResidualsBySaturations(
      ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxFunctionDerivatives,
      ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
      ConstMatrixRef pressureDerivativesX, ConstMatrixRef pressureDerivativesY,
      ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater,
      const Real timestep, const Real meshWidth
);


SparseMatrix deriveSaturationResidualsByPressures(const SparseMatrix& pressureSystem,
                                                               ConstMatrixRef fluxFunctionFactors,
                                                               ConstMatrixRef darcyVelocitiesX,
                                                               ConstMatrixRef darcyVelocitiesY,
                                                               ConstMatrixRef mobilities,
                                                               const Real timestep, const Real meshWidth);
SparseMatrix derivePressureResidualsByLogPermeabilities(ConstMatrixRef pressures, ConstMatrixRef totalMobilities);
SparseMatrix deriveSaturationResidualsByLogPermeabilities(ConstMatrixRef pressureGradientsX,
                                                                      ConstMatrixRef pressureGradientsY,
                                                                      ConstMatrixRef darcyVelocitiesX,
                                                                      ConstMatrixRef darcyVelocitiesY,
                                                                      ConstMatrixRef mobilities,
                                                                      ConstMatrixRef fluxFunctionFactors,
                                                                      const Real timestep, const Real meshWidth);

inline static Real hmeanDerivedBySecond(const Real a, const Real b) {
    return 2*std::pow(a/(a+b), 2);
}