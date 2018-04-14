//
// Created by Stefano Weidmann on 24.03.18.
//

#pragma once

#include "typedefs.hpp"
#include "cellindex.hpp"
#include "vectorToBeMappedAsMatrix.hpp"
#include "fixedParameters.hpp"
#include "simulationState.hpp"
#include "specialCells.hpp"



Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell);


SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities);

Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstVectorRef rhs);

void adaptRhsForPressure(const Real sourceAtDrillNow, VectorRef rhs, const int numberOfRows,
                         const int numberOfCols);

[[deprecated("Not used anymore by current formulation")]]
Vector projectSourcesIntoRange(ConstMatrixRef sources);

// Gradient zero at boundaries
Matrix computeXDerivative(ConstMatrixRef field, const Real meshWidth);

Matrix computeYDerivative(ConstMatrixRef field, const Real meshWidth);


void adaptPressureGradientsAtWell(const Real inflowNow, ConstMatrixRef mobilities, ConstMatrixRef pressures, MatrixRef pressureDerivativesX, MatrixRef pressureDerivativesY, const Real meshWidth);

Real getDerivativeAtCellBorder(CellIndex cell,
                               ConstMatrixRef xDerivative, ConstMatrixRef yDerivative,
                               const CellIndex::Direction whichBorder);
CellIndex centerIndexToBorderIndex(CellIndex centerIndex, const CellIndex::Direction whichBorder);
CellIndex borderIndexToCenterIndex(CellIndex borderIndex, const CellIndex::Direction whichBorder);

CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const int numberOfRows,
      const int numberOfCols);

Matrix computeTotalDarcyVelocitiesX(ConstMatrixRef totalMobilities, Matrix pressureDerivativesX);

Matrix computeTotalDarcyVelocitiesY(ConstMatrixRef totalMobilities, Matrix pressureDerivativesY);

Matrix computeFluxFunctionFactors(ConstMatrixRef saturations, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);

Matrix computeFluxesX(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesX);

Matrix computeFluxesY(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesY);


Matrix computeSaturationDivergences(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, const Real meshWidth);

MinimizationState doAMinimizerStep(MinimizationState oldState, ConstMatrixRef sensitivity);


Real clamp(const Real x, const Real minVal, const Real maxVal);
Matrix clamp(ConstMatrixRef x, const Real minVal, const Real maxVal);

Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater);

void advanceSaturationsInTime(const FixedParameters& params, MatrixRef saturationsWater,
                              ConstMatrixRef pressures, ConstMatrixRef totalMobilities, Real& time);

void stepForwardProblem(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& currentState, VectorRef pressureRhs);