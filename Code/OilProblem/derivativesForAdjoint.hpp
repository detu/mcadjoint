#pragma once
#include "typedefs.hpp"

SparseMatrix computePressureResidualsDerivedByPressure(const SparseMatrix& pressureSystem);
Matrix computeTotalMobilitiesDerivedBySaturationsWater(ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const Real dynamicViscosityOil, const Real dynamicViscosityWater);
SparseMatrix computePressureResidualsDerivedBySaturationWater(ConstMatrixRef pressures, ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater);
Matrix computeFluxFunctionFactorDerivatives(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);
SparseMatrix computeSaturationWaterResidualsDerivedBySaturationWater(ConstMatrixRef fluxFunctionDerivatives,
                                                                     ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
                                                                     const Real timestep, const Real meshWidth
);
SparseVector computeCostFunctionDerivedByPressure(const Real computedPressureAtDrill, const Real measuredPressureAtDrill, const int numberOfRows, const int numberOfCols);
SparseVector computeCostFunctionDerivedBySaturationsWater(const int numberOfRows, const int numberOfCols);

SparseVector computePressurePartOfDiagonalBlockTimesCostDerivedByState(const SparseMatrix& pressureResidualsByPressure, ConstMatrixRef pressures, const Real measuredPressureAtDrill);
SparseMatrix computeSaturationWaterResidualsDerivedByPressure(const SparseMatrix& pressureSystem, ConstMatrixRef fluxFunctionFactors,
                                                              ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
                                                              ConstMatrixRef mobilities,
                                                              const Real timestep, const Real meshWidth);
SparseMatrix computePressureResidualsByLogPermeability(ConstMatrixRef pressures, ConstMatrixRef totalMobilities);
SparseMatrix computeSaturationsWaterResidualsByLogPermeability(ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, ConstMatrixRef mobilities, const Real timestep, const Real meshWidth);

