//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_OILPROBLEM_HPP
#define STEFCOMMONHEADERS_OILPROBLEM_HPP

#include "typedefs.hpp"
#include "cellindex.hpp"

#if !defined(__GNUC__) && !defined(__attribute__)
    #define __attribute__(ignored)
#endif

__attribute__((pure))
Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell);


__attribute__((pure))
SparseMatrix assembleTransmissibilityMatrix(ConstMatrixRef totalMobilities);


__attribute__((pure))
VectorToBeMappedAsMatrix solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstMatrixRef sources);

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
      const int numberOfRows);

#endif //STEFCOMMONHEADERS_OILPROBLEM_HPP
