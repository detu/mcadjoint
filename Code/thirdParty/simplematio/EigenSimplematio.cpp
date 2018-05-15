#if EIGEN_WORLD_VERSION != 3
#warning "Only Eigen 3 supported!"
#endif

#include <stdexcept>
#include <cassert>
#include <string>
#include <fstream>
#include <type_traits>
#include <cstddef>


/**************************************************************************************************
 *                                  CONSTRUCTORS AND DESTRUCTOR                                   *
 **************************************************************************************************/

// mode_M is one of:
// - 'a': append (write, but do not delete previous contents)
// - 'w': write (delete previous contents)
// - 'r': read

SMIO::EigenMatFile::EigenMatFile(
    const char* filename
):
    mode_M('\0'),
    filename_M(nullptr),
    matFile_M(nullptr) {

    SMIO_SetErrorHandler(SMIO::EigenMatFile::rethrowSMIOError);

    open(filename);
}

SMIO::EigenMatFile::EigenMatFile(
):
    mode_M('\0'),
    filename_M(nullptr),
    matFile_M(nullptr) {
    SMIO_SetErrorHandler(SMIO::EigenMatFile::rethrowSMIOError);
}

SMIO::EigenMatFile::~EigenMatFile(
) {
    close();
}


/**************************************************************************************************
 *                                 EXPLICITLY OPENING AND CLOSING                                 *
 **************************************************************************************************/

void
SMIO::EigenMatFile::open(
    const char* filename
) {
    if (matFile_M != nullptr) {
        close();
    }
    filename_M = filename;
    if (std::ifstream(filename).good()) {
        mode_M = 'r';
    } else if (std::ofstream(filename).good()) {
        mode_M = 'w';
    } else {
        throw std::runtime_error(
            "Mat file \"" +
            std::string(filename) +
            "\" is neither readable nor writable!"
        );
    }

    matFile_M = SMIO_Open(filename, mode_M);
    if (mode_M == 'r') {
        populateAvailableVariables();
    }
}

void
SMIO::EigenMatFile::close(
) {
    if (matFile_M != nullptr) {
        SMIO_Close(matFile_M);
        matFile_M = nullptr;
    }
}

/**************************************************************************************************
 *                                       DIMENSION CHECKING                                       *
 **************************************************************************************************/

template <typename EigenMatrixType>
void
SMIO::EigenMatFile::checkDimensions(
    const std::size_t rowsAtRunTime,
    const std::size_t colsAtRunTime
) {

    constexpr std::size_t rowsAtCompileTime(Eigen::internal::traits<EigenMatrixType>::RowsAtCompileTime);
    constexpr std::size_t colsAtCompileTime(Eigen::internal::traits<EigenMatrixType>::RowsAtCompileTime);

    if ((rowsAtCompileTime != Eigen::Dynamic &&rowsAtRunTime != rowsAtCompileTime) || (colsAtCompileTime != Eigen::Dynamic && colsAtCompileTime != colsAtRunTime)) {
        throw std::runtime_error(
            "Expected to read a " + std::to_string(rowsAtCompileTime) + " x " + std::to_string(colsAtCompileTime)
            + " matrix but found a " + std::to_string(rowsAtRunTime) + " x " + std::to_string(colsAtRunTime)
        );
    }
}


/**************************************************************************************************
 *                                            READING                                             *
 **************************************************************************************************/

// read a variable from a mat file;
// reopens the file in read mode_M if necessary

template <typename EigenMatrixType>
EigenMatrixType
SMIO::EigenMatFile::readVariable(
    const char* variableName
) {
    assertReadyForReading(variableName);
    return internalReadVariable<EigenMatrixType>(variableName);
};

// dense and real
template <typename EigenMatrixType>
std::enable_if_t<
    SMIO::isRealMatrix<EigenMatrixType>::value && !SMIO::isSparseMatrix<EigenMatrixType>::value,
    EigenMatrixType
>
SMIO::EigenMatFile::internalReadVariable(
    const char* variableName
) {
    SMIO_DenseRealMatrix_t readMatrix(SMIO_ReadDenseRealMatrix(matFile_M, variableName));
    checkDimensions<EigenMatrixType>(readMatrix.rows, readMatrix.cols);
    return EigenMatrixType(
               std::move(Eigen::Map<SMIO::Matrix>(readMatrix.data, readMatrix.rows, readMatrix.cols))
           );
}

// dense and complex
template <typename EigenMatrixType>
std::enable_if_t<
    !SMIO::isRealMatrix<EigenMatrixType>::value && !SMIO::isSparseMatrix<EigenMatrixType>::value,
    EigenMatrixType
>
SMIO::EigenMatFile::internalReadVariable(
    const char* variableName
) {
    SMIO_DenseComplexMatrix_t readMatrix(SMIO_ReadDenseComplexMatrix(matFile_M, variableName));
    checkDimensions<EigenMatrixType>(readMatrix.rows, readMatrix.cols);
    return EigenMatrixType(
               std::move(Eigen::Map<SMIO::ComplexMatrix>(reinterpret_cast<SMIO::ComplexScalar*>(readMatrix.data), readMatrix.rows, readMatrix.cols))
           );
}

// sparse and real
template <typename EigenMatrixType>
std::enable_if_t<
    SMIO::isRealMatrix<EigenMatrixType>::value && SMIO::isSparseMatrix<EigenMatrixType>::value,
    EigenMatrixType
>
SMIO::EigenMatFile::internalReadVariable(
    const char* variableName
) {
    SMIO_SparseRealMatrix_t matrix(SMIO_ReadSparseRealMatrix(matFile_M, variableName));
    checkDimensions<EigenMatrixType>(matrix.rows, matrix.cols);

    #if EIGEN_MAJOR_VERSION >= 3
    return EigenMatrixType(
               std::move(
                   Eigen::Map<SMIO::SparseMatrix>(
                       matrix.rows,
                       matrix.cols,
                       matrix.numberOfNonzeroEntries,
                       matrix.colStarts,
                       matrix.rowIndices,
                       matrix.data,
                       nullptr
                   )
               )
           );
    #else
    return EigenMatrixType(
               std::move(
                   Eigen::MappedSparseMatrix<SMIO_RealScalar_t, SMIO::SparseMatrix::Flags>(
                       matrix.rows,
                       matrix.cols,
                       matrix.numberOfNonzeroEntries,
                       matrix.colStarts,
                       matrix.rowIndices,
                       matrix.data
                   )
               )
           );
    #endif
}

// sparse and complex
template <typename EigenMatrixType>
std::enable_if_t<
    !SMIO::isRealMatrix<EigenMatrixType>::value && SMIO::isSparseMatrix<EigenMatrixType>::value,
    EigenMatrixType
>
SMIO::EigenMatFile::internalReadVariable(
    const char* variableName
) {
    SMIO_SparseComplexMatrix_t matrix(SMIO_ReadSparseComplexMatrix(matFile_M, variableName));
    checkDimensions<EigenMatrixType>(matrix.rows, matrix.cols);


    #if EIGEN_MAJOR_VERSION >= 3
    return EigenMatrixType(
               std::move(
                   Eigen::Map<SMIO::SparseComplexMatrix>(
                       matrix.rows,
                       matrix.cols,
                       matrix.numberOfNonzeroEntries,
                       matrix.colStarts,
                       matrix.rowIndices,
                       reinterpret_cast<SMIO::ComplexScalar*>(matrix.data),
                       nullptr
                   )
               )
           );
    #else
    return EigenMatrixType(
               std::move(
                   Eigen::MappedSparseMatrix<SMIO::ComplexScalar, SMIO::SparseMatrix::Flags>(
                       matrix.rows,
                       matrix.cols,
                       matrix.numberOfNonzeroEntries,
                       matrix.colStarts,
                       matrix.rowIndices,
                       reinterpret_cast<SMIO::ComplexScalar*>(matrix.data)
                   )
               )
           );
    #endif
}

// scalar
template <typename EigenMatrixType>
std::enable_if_t<
    std::is_same<EigenMatrixType, SMIO::Scalar>::value || std::is_same<EigenMatrixType, SMIO::ComplexScalar>::value,
    EigenMatrixType
>
SMIO::EigenMatFile::internalReadVariable(
    const char* variableName
) {
    using Matrix1 = Eigen::Matrix<EigenMatrixType, 1, 1>;
    const Matrix1 scalar(internalReadVariable<Matrix1>(variableName));
    return scalar(0, 0);
}


/**************************************************************************************************
 *                                            WRITING                                             *
 **************************************************************************************************/

// truncate (delete contents) MatFile and open it for writing
void
SMIO::EigenMatFile::truncate(
) {
    availableVariables_M.clear();
    reopen('w');
}

// write a variable to a mat file;
void
SMIO::EigenMatFile::writeVariable(
    const char* variableName,
    SMIO::ConstMatrixRef matrix
) {
    assertReadyForWriting(variableName);
    deleteIfAlreadyPresent(variableName);

    SMIO_DenseRealMatrix_t matrixToWrite = SMIO_MakeDenseRealMatrix(
            matrix.rows(),
            matrix.cols(),
            const_cast<SMIO::Scalar*>(matrix.data())
                                           );

    SMIO_WriteDenseRealMatrix(
        matFile_M, variableName, &matrixToWrite
    );
}

void
SMIO::EigenMatFile::writeVariable(
    const char* variableName,
    SMIO::ConstComplexMatrixRef matrix
) {
    assertReadyForWriting(variableName);

    SMIO_DenseComplexMatrix_t matrixToWrite = SMIO_MakeDenseComplexMatrix(
                matrix.rows(),
                matrix.cols(),
                reinterpret_cast<SMIO_ComplexScalar_t*>(const_cast<SMIO::ComplexScalar*>(matrix.data()))
            );

    SMIO_WriteDenseComplexMatrix(
        matFile_M, variableName, &matrixToWrite
    );
}

void
SMIO::EigenMatFile::writeVariable(
    const char* variableName,
    SMIO::ConstSparseMatrixRef matrix
) {
    assertReadyForWriting(variableName);
    assert(matrix.isCompressed());

    const SMIO_SparseRealMatrix_t matrixToWrite {
        size_t(matrix.rows()),
        size_t(matrix.cols()),
        const_cast<SMIO_RealScalar_t*>(matrix.valuePtr()),
        size_t(matrix.nonZeros()),
        const_cast<int*>(matrix.innerIndexPtr()),
        int(matrix.nonZeros()),
        const_cast<int*>(matrix.outerIndexPtr()),
        int(matrix.outerSize() + 1)
    };

    SMIO_WriteSparseRealMatrix(
        matFile_M, variableName, &matrixToWrite
    );
}

void
SMIO::EigenMatFile::writeVariable(
    const char* variableName,
    SMIO::ConstSparseComplexMatrixRef matrix
) {
    assertReadyForWriting(variableName);
    assert(matrix.isCompressed());
    const SMIO_SparseComplexMatrix_t matrixToWrite {
        size_t(matrix.rows()),
        size_t(matrix.cols()),
        reinterpret_cast<SMIO_ComplexScalar_t*>(const_cast<SMIO::ComplexScalar*>((matrix.valuePtr()))),
        size_t(matrix.nonZeros()),
        const_cast<int*>(matrix.innerIndexPtr()),
        int(matrix.nonZeros()),
        const_cast<int*>(matrix.outerIndexPtr()),
        int(matrix.outerSize() + 1)
    };

    SMIO_WriteSparseComplexMatrix(
        matFile_M, variableName, &matrixToWrite
    );
}


void
SMIO::EigenMatFile::writeVariable(
    const char* variableName,
    const SMIO::Scalar value
) {
    return writeVariable(variableName, Eigen::Matrix<SMIO::Scalar, 1, 1>::Constant(value));
}

void
SMIO::EigenMatFile::writeVariable(
    const char* variableName,
    const SMIO::ComplexScalar value
) {
    return writeVariable(variableName, Eigen::Matrix<SMIO::ComplexScalar, 1, 1>::Constant(value));
}





/**************************************************************************************************
 *                                            INTERNAL                                            *
 **************************************************************************************************/

void
SMIO::EigenMatFile::populateAvailableVariables(
) {
    availableVariables_M.clear();
    Mat_Rewind(matFile_M);
    matvar_t* currentMatVarPointer;
    for (;;) {
        currentMatVarPointer = Mat_VarReadNextInfo(matFile_M);
        if (currentMatVarPointer == nullptr) {
            return;
        }

        if (!availableVariables_M.emplace(currentMatVarPointer->name).second) {
            throw std::runtime_error(
                "There is more than one variable of the name \"" +
                std::string(currentMatVarPointer->name) +
                "\" in the mat file \"" +
                std::string(filename_M) +
                "\"! Can't handle this!"
            );
        }
        Mat_VarFree(currentMatVarPointer);
    }
    Mat_Rewind(matFile_M);
}


bool
SMIO::EigenMatFile::deleteIfAlreadyPresent(
    const std::string& variableName
) {
    if (availableVariables_M.count(variableName) == 1) {
        // "deletes" by copying all other variables into a new file --> expensive

        Mat_VarDelete(matFile_M, variableName.c_str());
        // FIXME: Why do i need to reopen the file?
        reopen('a');
        return true;
    }
    return false;
}

void
SMIO::EigenMatFile::reopen(
    const char newmode_M
) {
    SMIO_Close(matFile_M);
    matFile_M = SMIO_Open(filename_M, newmode_M);
    mode_M = newmode_M;
}


void
SMIO::EigenMatFile::assertReadyForReading(
    const std::string& variableName
) {
    if (matFile_M == nullptr) {
        throw std::logic_error("No open file! Open a file with .open() before reading!");
    } else if (mode_M != 'r') {
        reopen('r');
    }

    if (availableVariables_M.count(variableName) == 0) {
        throw std::runtime_error(
            "Wanted to read variable \"" +
            variableName +
            "\" from file \"" +
            std::string(filename_M) +
            " but there's no such variable there!"
        );
    }
}

void
SMIO::EigenMatFile::assertReadyForWriting(
    const std::string& variableName
) {
    if (matFile_M == nullptr) {
        throw std::logic_error("No open file! Open a file with .open() before writing!");
    } else if (mode_M != 'w' && mode_M != 'a') {
        reopen('a');
    }


    if (!deleteIfAlreadyPresent(variableName)) {
        // variable name is new
        // add it to the list of available variables
        availableVariables_M.insert(variableName);
    }

}


void
SMIO::EigenMatFile::rethrowSMIOError(
    const char* errorMessage
) {
    throw std::runtime_error(errorMessage);
}

