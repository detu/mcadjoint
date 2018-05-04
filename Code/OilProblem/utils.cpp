//
// Created by Stefano Weidmann on 13.03.18.
//

#include <algorithm>
#include "utils.hpp"

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