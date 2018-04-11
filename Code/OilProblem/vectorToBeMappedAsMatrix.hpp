//
// Created by Stefano Weidmann on 02.04.18.
//

#pragma once
#include "typedefs.hpp"
struct VectorToBeMappedAsMatrix {
    Vector vec;
    Eigen::Map<Matrix> map;

    inline VectorToBeMappedAsMatrix(Vector input, const int matrixRows, const int matrixCols):
          vec(std::move(input)), map(vec.data(), matrixRows, matrixCols) {}

    inline VectorToBeMappedAsMatrix(const VectorToBeMappedAsMatrix& other):
          VectorToBeMappedAsMatrix(other.vec, other.map.rows(), other.map.cols()) {}

    inline VectorToBeMappedAsMatrix(VectorToBeMappedAsMatrix&& other):
          VectorToBeMappedAsMatrix(std::move(other.vec), other.map.rows(), other.map.cols()) {}

    inline VectorToBeMappedAsMatrix(const int matrixRows, const int matrixCols):
          vec(matrixRows * matrixCols), map(vec.data(), matrixRows, matrixCols) {}


    inline VectorToBeMappedAsMatrix& operator =(ConstVectorRef vector) {
        vec = vector;
        map = Eigen::Map<Matrix>(vec.data(), map.rows(), map.cols());
        return *this;
    }
};


struct MinimizationState {
    Real stepSize;
    Real cost;
    VectorToBeMappedAsMatrix parameters;

    inline MinimizationState(const int matrixRows, const int matrixCols):
          stepSize(NAN), cost(NAN), parameters(matrixRows, matrixCols) {}
};


