//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_TYPEDEFS_HPP
#define STEFCOMMONHEADERS_TYPEDEFS_HPP
#include <Eigen/Sparse>
#include <Eigen/Dense>

using Real = double;
using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using VectorRef = Eigen::Ref<Vector>;
using ConstVectorRef = const Eigen::Ref<const Vector>;

using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixRef = Eigen::Ref<Matrix>;
using ConstMatrixRef = const Eigen::Ref<const Matrix>;

using SparseMatrix = Eigen::SparseMatrix<Real, Eigen::ColMajor>;



#endif //STEFCOMMONHEADERS_TYPEDEFS_HPP