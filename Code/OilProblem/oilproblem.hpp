//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_OILPROBLEM_HPP
#define STEFCOMMONHEADERS_OILPROBLEM_HPP

#include "typedefs.hpp"
#include "cellindex.hpp"
Real computeTransmissibility(ConstMatrixRef lambdas, const CellIndex& fromCell, const CellIndex& toCell);

CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const long numberOfRows);

SparseMatrix assembleTransmissibilityMatrix(ConstMatrixRef lambdas);

Matrix solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstMatrixRef sources);

#endif //STEFCOMMONHEADERS_OILPROBLEM_HPP
