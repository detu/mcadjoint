//
// Created by Stefano Weidmann on 30.04.18.
//
#include "dumpToMatFile.hpp"
#include <unordered_map>
#include <cstdlib>
#include <stefCommonHeaders/omp_mutex.hpp>
#include "logging.hpp"

using MatricesToDump = std::unordered_map<std::string, Matrix>;
using SparseMatricesToDump = std::unordered_map<std::string, SparseMatrix>;

static MatricesToDump MATRICES_TO_DUMP;
static SparseMatricesToDump SPARSE_MATRICES_TO_DUMP;

static std::string MAT_FILE_NAME = "";


static void dumpMatricesWithoutLock() {
    if (!MAT_FILE_NAME.empty()) {
        SMIO::EigenMatFile matFile(MAT_FILE_NAME.c_str());
        for (const std::pair<std::string, Matrix>& entry: MATRICES_TO_DUMP) {
            matFile.writeVariable(entry.first.c_str(), entry.second);
        }

        for (const auto& entry: SPARSE_MATRICES_TO_DUMP) {
            matFile.writeVariable(entry.first.c_str(), entry.second);
        }

        matFile.close();
    }
}

void writeToMatFile() {
    dumpMatricesWithoutLock();
}

void dumpInThisMatFile(const std::string& matFileName) {
    MAT_FILE_NAME = matFileName;
}






void dumpThis(const char* varName, ConstMatrixRef matrix) {
    Matrix& matrixToDump = MATRICES_TO_DUMP[std::string(varName)];
    matrixToDump.resizeLike(matrix);
    matrixToDump = matrix;


}

void dumpThis(const char* varName, const Real scalar) {
    return dumpThis(varName, Matrix::Constant(1, 1, scalar));
}

void dumpThis(const char* varName, const SparseMatrix& matrix) {
    SparseMatrix& matrixToDump = SPARSE_MATRICES_TO_DUMP[std::string(varName)];
    matrixToDump.resize(matrix.rows(), matrix.cols());
    matrixToDump.resizeNonZeros(matrix.nonZeros());

    matrixToDump = matrix;
}

void dumpThis(const char* varName, const std::vector<Real>& vector) {
    return dumpThis(varName, Eigen::Map<const Vector>(vector.data(), vector.size()));
}

void dumpThis(const char* varName, const std::vector<Matrix>& matrices) {
    if (matrices.empty()) {
        return;
    }

    const int rowsOfFirstMatrix = matrices[0].rows();

    int totalCols = 0;

    for (const Matrix& matrix: matrices) {
        const int rowsOfThisMatrix = matrix.rows();
        if (rowsOfThisMatrix != rowsOfFirstMatrix) {
            throw std::logic_error("Rows aren't constant!");
        }

        const int colsOfThisMatrix = matrix.cols();
        totalCols += colsOfThisMatrix;
    }

    Matrix megaMatrix(rowsOfFirstMatrix, totalCols);

    int beginColIndex = 0;

    for (const Matrix& matrix: matrices) {
        const int colsOfThisMatrix = matrix.cols();

        megaMatrix.middleCols(beginColIndex, colsOfThisMatrix) = matrix;

        beginColIndex += colsOfThisMatrix;
    }

    return dumpThis(varName, megaMatrix);
}