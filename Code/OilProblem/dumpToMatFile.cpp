//
// Created by Stefano Weidmann on 30.04.18.
//
#include "dumpToMatFile.hpp"
#include <unordered_map>
#include <cstdlib>
#include <mutex>
#include <csignal>

using MatricesToDump = std::unordered_map<std::string, Matrix>;
using SparseMatricesToDump = std::unordered_map<std::string, SparseMatrix>;

static MatricesToDump MATRICES_TO_DUMP;
static SparseMatricesToDump SPARSE_MATRICES_TO_DUMP;

static std::string MAT_FILE_NAME = "";


static std::mutex MAT_MUTEX;


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
    std::lock_guard<std::mutex> lockGuard(MAT_MUTEX);
    dumpMatricesWithoutLock();
}

void dumpInThisMatFile(const std::string& matFileName) {
    std::lock_guard<std::mutex> lockGuard(MAT_MUTEX);
    MAT_FILE_NAME = matFileName;

}






void dumpThis(const char* varName, ConstMatrixRef matrix) {
    std::lock_guard<std::mutex> lockGuard(MAT_MUTEX);

    Matrix& matrixToDump = MATRICES_TO_DUMP[std::string(varName)];
    matrixToDump.resizeLike(matrix);
    matrixToDump = matrix;


}

void dumpThis(const char* varName, const Real scalar) {
    return dumpThis(varName, Matrix::Constant(1, 1, scalar));
}

void dumpThis(const char* varName, const SparseMatrix& matrix) {
    std::lock_guard<std::mutex> lockGuard(MAT_MUTEX);

    SparseMatrix& matrixToDump = SPARSE_MATRICES_TO_DUMP[std::string(varName)];
    matrixToDump.resize(matrix.rows(), matrix.cols());
    matrixToDump.resizeNonZeros(matrix.nonZeros());

    matrixToDump = matrix;
}