//
// Created by Stefano Weidmann on 24.03.18.
//

#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/OrderingMethods>
#include <stefCommonHeaders/xoroshiro.h>
#include <functional>

using Real = double;

using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using VectorRef = Eigen::Ref<Vector>;
using ConstVectorRef = const Eigen::Ref<const Vector>&;

using SignalHandler = void (*)(int);

using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixRef = Eigen::Ref<Matrix>;
using ConstMatrixRef = const Eigen::Ref<const Matrix>&;

using CostFunction = std::function<Real(ConstMatrixRef)>;

using WellFunction = std::function<Real(Real)>;
using PressureFunction = WellFunction;

using Rng = XoroshiroRandomNumberEngine;

struct Point {
    Real x;
    Real y;
};


using SparseMatrix = Eigen::SparseMatrix<Real, Eigen::ColMajor>;


