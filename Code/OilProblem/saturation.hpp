#pragma once
#include "typedefs.hpp"
#include "fixedParameters.hpp"

bool advanceSaturationsInTime(const FixedParameters& params, MatrixRef saturationsWater, ConstMatrixRef pressures, ConstMatrixRef totalMobilities, Real& time);
Real computeTimestep(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY, const Real meshWidth, const Real finalTime, const Real time);
Real getFirstTimestep();
Matrix computeSaturationDivergences(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, const Real meshWidth);
Real clamp(const Real x, const Real minVal, const Real maxVal);
Matrix clamp(ConstMatrixRef x, const Real minVal, const Real maxVal);
Matrix computeFluxesY(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesY);
Matrix computeFluxesX(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesX);
Matrix computeFluxFunctionFactors(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);
