//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_OILPROBLEM_HPP
#define STEFCOMMONHEADERS_OILPROBLEM_HPP

#include "typedefs.hpp"
#include "cellindex.hpp"
#include "vectorToBeMappedAsMatrix.hpp"

#if !defined(__GNUC__) && !defined(__attribute__)
    #define __attribute__(ignored)
#endif

__attribute__((pure))
Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell);


__attribute__((pure))
SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities);

__attribute__((pure))
Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstVectorRef negatedSourcesProjectedIntoRange, ConstVectorRef pressureGuess);

__attribute__((pure))
Vector augmentSources(ConstMatrixRef sources);

__attribute__((pure))
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
Matrix computeTotalDarcyVelocitiesX(ConstMatrixRef totalTransmissibilities, Matrix pressureDerivativesX);

__attribute__((pure))
Matrix computeTotalDarcyVelocitiesY(ConstMatrixRef totalTransmissibilities, Matrix pressureDerivativesY);

__attribute__((pure))
Matrix computeFluxFunctionFactors(ConstMatrixRef saturations, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil);

__attribute__((pure))
Matrix computeFluxesX(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesX);

__attribute__((pure))
Matrix computeFluxesY(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesY);

__attribute__((pure))
Matrix advanceStateInTime(Matrix state, ConstMatrixRef derivativeInTime);

__attribute__((pure))
Matrix computeSaturationDivergence(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, const Real meshWidth);

__attribute__((pure))
MinimizationState doAMinimizerStep(MinimizationState oldState, ConstMatrixRef sensitivity);

__attribute__((pure))
Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater);


#endif //STEFCOMMONHEADERS_OILPROBLEM_HPP
