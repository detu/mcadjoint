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


bool allOf(const SparseMatrix& matrix, const Predicate<Real>& predicate);
bool allOf(ConstMatrixRef matrix, const Predicate<Real>& predicate);
bool anyOf(const SparseMatrix& matrix, const Predicate<Real>& predicate);
bool anyOf(ConstMatrixRef matrix, const Predicate<Real>& predicate);

bool allFinite(ConstMatrixRef matrix);
bool allFinite(const SparseMatrix& matrix);
