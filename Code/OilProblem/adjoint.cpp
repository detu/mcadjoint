//
// Created by Stefano Weidmann on 17.04.18.
//

#include "oilproblem.hpp"

SparseMatrix makeMatrixColStochastic(SparseMatrix matrix) {
    for (int col = 0; col < matrix.cols(); ++col) {
        matrix.col(col) /= matrix.col(col).sum();
    }

    return matrix;
}

Real getTransitionProbability(const int fromState, const int toState, const SparseMatrix& transitionProbabilitiesColStochastic) {
    return transitionProbabilitiesColStochastic.coeff(toState, fromState);
}

