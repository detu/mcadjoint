//
// Created by Stefano Weidmann on 13.03.18.
//

#pragma once
#include <string>
#include "typedefs.hpp"


std::string lowerCase(std::string str);
Real frobeniusNormSquared(const SparseMatrix& matrix);
DiagonalMatrix extractInverseDiagonalMatrix(const SparseMatrix& matrix);
Real sumOfAbsEntries(const SparseMatrix& matrix);

static inline int sign(const Real number) {
    return 2 * (number > 0) - 1;
}
