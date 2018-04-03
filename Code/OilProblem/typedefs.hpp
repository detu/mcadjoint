//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_TYPEDEFS_HPP
#define STEFCOMMONHEADERS_TYPEDEFS_HPP
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <functional>

using Real = double;
using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using VectorRef = Eigen::Ref<Vector>;
using ConstVectorRef = const Eigen::Ref<const Vector>&;

using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixRef = Eigen::Ref<Matrix>;
using ConstMatrixRef = const Eigen::Ref<const Matrix>&;

using SparseVector = Eigen::SparseVector<Real>;

using CostFunction = std::function<Real(Matrix)>;

using WellFunction = std::function<Real(Real)>;
using PressureFunction = WellFunction;



struct Point {
    Real x;
    Real y;
};


#ifdef USE_PARDISO
#include <Eigen/PardisoSupport>
using SparseMatrix = Eigen::SparseMatrix<Real, Eigen::RowMajor>;
using SparseSolver = Eigen::PardisoLU<SparseMatrix>;
#else
using SparseMatrix = Eigen::SparseMatrix<Real, Eigen::ColMajor>;
using SparseSolver = Eigen::SparseLU<SparseMatrix>;

#endif


#endif //STEFCOMMONHEADERS_TYPEDEFS_HPP
