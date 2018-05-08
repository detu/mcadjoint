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

Real computeCostFunctionDerivedByPressureEntry(const Real computedPressureAtDrill, const Real measuredPressureAtDrill,
                                               const int numberOfRows, const int numberOfCols,
                                               const CellIndex& derivedByCell);
Real computeCostFunctionDerivedBySaturationsWaterEntry(const int numberOfRows, const int numberOfCols,
                                                       const CellIndex& derivedByCell);

SparseMatrix computeSaturationWaterResidualsDerivedByPressure(const SparseMatrix& pressureSystem, ConstMatrixRef fluxFunctionFactors,
                                                              ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
                                                              ConstMatrixRef mobilities,
                                                              const Real timestep, const Real meshWidth);
SparseMatrix computePressureResidualsByLogPermeability(ConstMatrixRef pressures, ConstMatrixRef totalMobilities);
SparseMatrix computeSaturationsWaterResidualsByLogPermeability(ConstMatrixRef pressureGradientsX, ConstMatrixRef pressureGradientsY, ConstMatrixRef mobilities, const Real timestep, const Real meshWidth);
