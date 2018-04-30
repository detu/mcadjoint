#pragma once
#include "typedefs.hpp"
#include "cellindex.hpp"

Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities,
                                   ConstVectorRef& rhs);
Vector computeRhsForPressureSystem(const Real sourceAtDrillNow, const int numberOfRows, const int numberOfCols);
SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities);
CellIndex pressureToTransmissibilityIndex(const CellIndex& fromCell, const CellIndex& toCell,const int numberOfRows);
Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell);
Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater);