//
// Created by Stefano Weidmann on 16.05.18.
//

#pragma once

#include "typedefs.hpp"

enum class WhichPreconditioner {
    NONE, DIAGONAL, PRESSURE_BY_PRESSURE, Q_TRANSPOSE_FROM_QR_AND_DIAGONAL, CHOLESKY_L_INV
};

inline static bool preconditionerLeavesJustOnesOnDiagonal(const WhichPreconditioner& whichPreconditioner) {
    return whichPreconditioner != WhichPreconditioner::NONE;
}

void preconditionMatrices(
      SparseMatrix& pressuresByPressures, SparseMatrix& saturationsByPressures,
      SparseMatrix& pressuresBySaturations, SparseMatrix& saturationsBySaturations,
      VectorRef b, const PressureSolver& pressureSolver, const WhichPreconditioner whichPreconditioner);