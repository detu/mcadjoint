//
// Created by Stefano Weidmann on 13.03.18.
//

#include <algorithm>
#include "utils.hpp"
#include "logging.hpp"

std::string lowerCase(std::string str) {
    transform(str.begin(), str.end(), str.begin(), [] (char c) -> char {return char(tolower(c));});
    return str;
}

Real frobeniusNormSquared(const SparseMatrix& matrix) {
    const Real* beginValues = matrix.valuePtr();
    const Real* endValues = beginValues + matrix.nonZeros();

    Real norm = 0;

    for (const Real* entryPtr = beginValues; entryPtr < endValues; ++entryPtr) {
        norm += std::pow(*entryPtr, 2);
    }

    return norm;
}

DiagonalMatrix extractInverseDiagonalMatrix(const SparseMatrix& matrix) {
    const DiagonalMatrix diagonalMatrix = matrix.diagonal().asDiagonal();

    constexpr bool printInverseDiagonal = false;

    if (printInverseDiagonal) {
        LOGGER->debug("inverse diagonal", Vector(diagonalMatrix.inverse()));
    }
    return diagonalMatrix.inverse();
}

Real sumOfAbsEntries(const SparseMatrix& matrix) {
    const Real* begin = matrix.valuePtr();
    const Real* end = begin + matrix.nonZeros();

    Real absSum = 0;
    for (const Real* entryPtr = begin; entryPtr < end; ++entryPtr) {
        absSum += std::abs(*entryPtr);
    }

    return absSum;
}