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
        vec.resizeLike(vector);
        vec = vector;
        map = Eigen::Map<Matrix>(vec.data(), map.rows(), map.cols());
        return *this;
    }

    /*inline operator Eigen::Ref<const Matrix>() const {
        return map;
    }

    inline operator Eigen::Ref<Matrix>() {
        return map;
    }

    inline operator Eigen::Ref<const Vector>() const {
        return vec;
    }

    inline operator Eigen::Ref<Vector>() {
        return vec;
    }*/

};





