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
#include "randomWalkState.hpp"
#include "bVectorSurrogate.hpp"
#include "cMatrixSurrogate.hpp"
#include "minimizationState.hpp"

Real getFirstTimestep();

Matrix computeFluxFunctionFactorDerivatives(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);

SparseMatrix computePressureResidualsDerivedByPressure(const SparseMatrix& pressureSystem);
SparseMatrix computePressureResidualsDerivedBySaturationWater(ConstMatrixRef pressures, ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater);
SparseMatrix computeSaturationWaterResidualsDerivedBySaturationWater(ConstMatrixRef fluxFunctionDerivatives,
                                                                     ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
                                                                     const Real timestep, const Real meshWidth
);
SparseMatrix computeSaturationWaterResidualsDerivedByPressure(const SparseMatrix& pressureSystem, ConstMatrixRef fluxFunctionFactors,
                                                              ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
                                                              ConstMatrixRef mobilities,
                                                              const Real timestep, const Real meshWidth);

Matrix computeTotalMobilitiesDerivedBySaturationsWater(ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const Real dynamicViscosityOil, const Real dynamicViscosityWater);

Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell);


SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities);

SparseMatrix computePressureResidualsByLogPermeability(ConstMatrixRef pressures, ConstMatrixRef totalMobilities);
SparseMatrix computeSaturationsWaterResidualsByLogPermeability(ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, ConstMatrixRef mobilities, const Real timestep, const Real meshWidth);

Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstVectorRef rhs);

void adaptRhsForPressure(const Real sourceAtDrillNow, VectorRef rhs, const int numberOfRows,
                         const int numberOfCols);

[[deprecated("Not used anymore by current formulation")]]
Vector projectSourcesIntoRange(ConstMatrixRef sources);

// Gradient zero at boundaries
Matrix computeXDerivative(ConstMatrixRef field, const Real meshWidth);

Matrix computeYDerivative(ConstMatrixRef field, const Real meshWidth);

bool transitionState(RandomWalkState& currentState, const BVectorSurrogate& b,
                     const SparseMatrix& pressureResidualsByPressures,
                     const SparseMatrix& pressureResidualsBySaturationsWater,
                     const SparseMatrix& saturationsWaterResidualsByPressure,
                     const SparseMatrix& saturationsWaterResidualsBySaturationsWater,
                     const int numberOfRows, const int numberOfCols, Rng& rng);

std::vector<RandomWalkState> initializeRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters, const BVectorSurrogate& b, const CMatrixSurrogate& c);


void adaptPressureGradientsAtWell(const Real inflowNow, ConstMatrixRef mobilities, ConstMatrixRef pressures, MatrixRef pressureDerivativesX, MatrixRef pressureDerivativesY, const Real meshWidth);

Real getDerivativeAtCellBorder(CellIndex cell,
                               ConstMatrixRef xDerivative, ConstMatrixRef yDerivative,
                               const CellIndex::Direction whichBorder);
CellIndex centerIndexToBorderIndex(CellIndex centerIndex, const CellIndex::Direction whichBorder);
CellIndex borderIndexToCenterIndex(CellIndex borderIndex, const CellIndex::Direction whichBorder);

CellIndex pressureToTransmissibilityIndex(const CellIndex& fromCell, const CellIndex& toCell, const int numberOfRows);

Matrix computeTotalDarcyVelocitiesX(ConstMatrixRef totalMobilities, Matrix pressureDerivativesX);

Matrix computeTotalDarcyVelocitiesY(ConstMatrixRef totalMobilities, Matrix pressureDerivativesY);

Matrix computeFluxFunctionFactors(ConstMatrixRef saturations, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);

Matrix computeFluxesX(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesX);

Matrix computeFluxesY(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesY);


Matrix computeSaturationDivergences(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, const Real meshWidth);

MinimizationState doAMinimizerStep(MinimizationState oldState, ConstMatrixRef sensitivity);

SparseMatrix makeMatrixColStochastic(SparseMatrix matrix);

Real clamp(const Real x, const Real minVal, const Real maxVal);
Matrix clamp(ConstMatrixRef x, const Real minVal, const Real maxVal);

Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater);

bool advanceSaturationsInTime(const FixedParameters& params, MatrixRef saturationsWater,
                              ConstMatrixRef pressures, ConstMatrixRef totalMobilities, Real& time);

bool stepForwardProblem(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& currentState, VectorRef pressureRhs);