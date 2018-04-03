//
// Created by Stefano Weidmann on 02.04.18.
//

#ifndef STEFCOMMONHEADERS_VECTORTOBEMAPPEDASMATRIX_HPP
#define STEFCOMMONHEADERS_VECTORTOBEMAPPEDASMATRIX_HPP
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
};


struct MinimizationState {
    Real stepSize;
    Real cost;
    VectorToBeMappedAsMatrix parameters;

    inline MinimizationState(const int matrixRows, const int matrixCols):
          stepSize(NAN), cost(NAN), parameters(matrixRows, matrixCols) {}
};

struct SimulationState {
    Matrix saturationsWater;
    VectorToBeMappedAsMatrix pressures;

    inline SimulationState(const int matrixRows, const int matrixCols):
          saturationsWater(matrixRows, matrixCols), pressures(matrixRows, matrixCols) {}
};

#endif //STEFCOMMONHEADERS_VECTORTOBEMAPPEDASMATRIX_HPP
