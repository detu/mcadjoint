//
// Created by Stefano Weidmann on 24.03.18.
//

#pragma once

#include "typedefs.hpp"
#include "cellindex.hpp"
#include "vectorToBeMappedAsMatrix.hpp"
#include "fixedParameters.hpp"
#include "simulationState.hpp"

#if !defined(__GNUC__) && !defined(__attribute__)
    #define __attribute__(ignored)
#endif

__attribute__((pure))
Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell);


__attribute__((pure))
SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities);

__attribute__((pure))
Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstVectorRef rhs, ConstVectorRef pressureGuess);

void adaptRhsForPressure(const Real sourceAtWellNow, const Real pressureAtDrillNow, VectorRef rhs, const int numberOfRows,
                         const int numberOfCols);

__attribute__((pure))
[[deprecated("Not used anymore by current formulation")]]
Vector projectSourcesIntoRange(ConstMatrixRef sources);

// Gradient zero at boundaries
__attribute__((pure))
Matrix computeXDerivative(ConstMatrixRef field, const Real meshWidth);

__attribute__((pure))
Matrix computeYDerivative(ConstMatrixRef field, const Real meshWidth);


__attribute__((pure))
Real getDerivativeAtCellBorder(CellIndex cell,
                               ConstMatrixRef xDerivative, ConstMatrixRef yDerivative,
                               const CellIndex::Direction whichBorder);

__attribute__((pure))
CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const int numberOfRows,
      const int numberOfCols);

__attribute__((pure))
Matrix computeTotalDarcyVelocitiesX(ConstMatrixRef totalMobilities, Matrix pressureDerivativesX);

__attribute__((pure))
Matrix computeTotalDarcyVelocitiesY(ConstMatrixRef totalMobilities, Matrix pressureDerivativesY);

__attribute__((pure))
Matrix computeFluxFunctionFactors(ConstMatrixRef saturations, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);

__attribute__((pure))
Matrix computeFluxesX(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesX);

__attribute__((pure))
Matrix computeFluxesY(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesY);


__attribute__((pure))
Matrix computeSaturationDivergences(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, const Real meshWidth);

__attribute__((pure))
MinimizationState doAMinimizerStep(MinimizationState oldState, ConstMatrixRef sensitivity);

__attribute__((pure))
Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater);

void advanceSaturationsInTime(const FixedParameters& params, MatrixRef saturationsWater,
                              ConstMatrixRef pressures, ConstMatrixRef totalMobilities, Real& time);

void stepForwardProblem(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& currentState, VectorRef pressureRhs);