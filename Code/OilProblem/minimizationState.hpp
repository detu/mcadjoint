//
// Created by stefanow on 4/27/18.
//

#pragma once
#include "vectorToBeMappedAsMatrix.hpp"
#include "typedefs.hpp"

struct MinimizationState {
    Real stepSize;
    Real cost;
    VectorToBeMappedAsMatrix parameters;

    inline MinimizationState(const int matrixRows, const int matrixCols):
          stepSize(NAN), cost(NAN), parameters(matrixRows, matrixCols) {}
};
