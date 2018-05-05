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
