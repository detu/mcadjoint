/*
88888888888 88   ,ad8888ba,  88888888888 888b      88
88          88  d8"'    `"8b 88          8888b     88
88          88 d8'           88          88 `8b    88
88aaaaa     88 88            88aaaaa     88  `8b   88
88"""""     88 88      88888 88"""""     88   `8b  88
88          88 Y8,        88 88          88    `8b 88
88          88  Y8a.    .a88 88          88     `8888
88888888888 88   `"Y88888P"  88888888888 88      `888

      a8"           "8a
    a8"               "8a
  a8"                   "8a
a8"                       "8a
"8a   aaaaaaaa aaaaaaaa   a8"
  "8a """""""" """""""" a8"
    "8a               a8"
      "8a           a8"

       88b           d88        db   888888888888
       888b         d888       d88b       88
       88`8b       d8'88      d8'`8b      88
       88 `8b     d8' 88     d8'  `8b     88
       88  `8b   d8'  88    d8YaaaaY8b    88
       88   `8b d8'   88   d8""""""""8b   88
888    88    `888'    88  d8'        `8b  88
888    88     `8'     88 d8'          `8b 88

*/

#pragma once
#include "simplematio.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_set>
#include <string>
#include <complex>
#include <cstddef>


#if __cplusplus < 201402L
namespace std {
    template <bool condition, typename T>
    using enable_if_t = typename std::enable_if<condition, T>::type;
}
#endif

namespace SMIO {

    
    /**************************************************************************************************
     **************************************************************************************************
     **************************************************************************************************
     ***                                          TYPEDEFS                                          ***
     **************************************************************************************************
     **************************************************************************************************
     **************************************************************************************************/

    using Scalar = SMIO_RealScalar_t;
    using ComplexScalar = std::complex<Scalar>;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>; // Eigen automatically transposes RowMajor matrices
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Eigen::ColMajor>; // Eigen automatically transposes RowMajor matrices
    using ComplexMatrix = Eigen::Matrix<ComplexScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using ComplexVector = Eigen::Matrix<ComplexScalar, Eigen::Dynamic, 1, Eigen::ColMajor>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor>;
    using SparseComplexMatrix = Eigen::SparseMatrix<ComplexScalar, Eigen::ColMajor>;
    using MatrixRef = Eigen::Ref<Matrix>;
    using ConstMatrixRef = const Eigen::Ref<const Matrix>&;
    using VectorRef = Eigen::Ref<Vector>;
    using ConstVectorRef = const Eigen::Ref<const Vector>&;
    using ComplexVectorRef = Eigen::Ref<ComplexVector>;
    using ConstComplexVectorRef = const Eigen::Ref<const ComplexVector>&;
    using ComplexMatrixRef = Eigen::Ref<ComplexMatrix>;
    using ConstComplexMatrixRef = const Eigen::Ref<const ComplexMatrix>&;

    // FIXME: Issue #1 (see https://bitbucket.org/HansiMeier/simplematio/issues/1/eigen-33-ambiguous-ref-for-sparse-matrices)
    #if EIGEN_MAJOR_VERSION >= 3 && defined(SIMPLEMATIO_ISSUE_ONE_RESOLVED)
    using SparseMatrixRef = Eigen::Ref<SparseMatrix, Eigen::StandardCompressedFormat>;
    using ConstSparseMatrixRef = const Eigen::Ref<const SparseMatrix, Eigen::StandardCompressedFormat>& ;
    using SparseComplexMatrixRef = Eigen::Ref<SparseComplexMatrix, Eigen::StandardCompressedFormat>;
    using ConstSparseComplexMatrixRef = const Eigen::Ref<const SparseComplexMatrix, Eigen::StandardCompressedFormat>& ;
    #else
    using SparseMatrixRef = SparseMatrix&;
    using ConstSparseMatrixRef = const SparseMatrix&; // before Eigen 3.3 there is no specialization of Eigen::Ref for sparse matrices
    using SparseComplexMatrixRef = SparseComplexMatrix&;
    using ConstSparseComplexMatrixRef = const SparseComplexMatrix&;
    #endif

    /**************************************************************************************************
     **************************************************************************************************
     **************************************************************************************************
     ***                                           TRAITS                                           ***
     **************************************************************************************************
     **************************************************************************************************
     **************************************************************************************************/

    /**************************************************************************************************
     *                                 CHECK IF MATRIX TYPE IS SPARSE                                 *
     **************************************************************************************************/

    template <typename EigenMatrixType>
    struct isSparseMatrix {
        const static bool value;
    };

    template <typename _Scalar, int _Options, typename _Index>
    struct isSparseMatrix<typename Eigen::SparseMatrix<_Scalar, _Options, _Index>> {
        const static bool value;
    };

    template <typename _Scalar, int _Options, typename _Index>
    const bool isSparseMatrix<typename Eigen::SparseMatrix<_Scalar, _Options, _Index>>::value(true);

    template <typename _Scalar, int _Rows, int _Cols, int _MaxRows, int _MaxCols>
    struct isSparseMatrix<Eigen::Matrix<_Scalar, _Rows, _Cols, Eigen::ColMajor, _MaxRows, _MaxCols>> {
        const static bool value;
    };

    template <typename _Scalar, int _Rows, int _Cols, int _MaxRows, int _MaxCols>
    const bool isSparseMatrix<Eigen::Matrix<_Scalar, _Rows, _Cols, Eigen::ColMajor, _MaxRows, _MaxCols>>::value(false);


    /**************************************************************************************************
     *                                  CHECK IF MATRIX TYPE IS REAL                                  *
     **************************************************************************************************/

    template <typename EigenMatrixType>
    struct isRealMatrix {
        const static bool value;
    };

    template <typename _Scalar, int _Options, typename _Index>
    struct isRealMatrix<typename Eigen::SparseMatrix<_Scalar, _Options, _Index>> {
        const static bool value;
    };

    template <typename _Scalar, int _Options, typename _Index>
    const bool isRealMatrix<typename Eigen::SparseMatrix<_Scalar, _Options, _Index>>::value(std::is_same<SMIO::Scalar, _Scalar>::value);

    template <typename _Scalar, int _Rows, int _Cols, int _MaxRows, int _MaxCols>
    struct isRealMatrix<Eigen::Matrix<_Scalar, _Rows, _Cols, Eigen::ColMajor, _MaxRows, _MaxCols>> {
        const static bool value;
    };

    template <typename _Scalar, int _Rows, int _Cols, int _MaxRows, int _MaxCols>
    const bool isRealMatrix<Eigen::Matrix<_Scalar, _Rows, _Cols, Eigen::ColMajor, _MaxRows, _MaxCols>>::value(std::is_same<SMIO::Scalar, _Scalar>::value);


    /**************************************************************************************************
     **************************************************************************************************
     **************************************************************************************************
     ***                                  EIGEN MAT FILE INTERFACE                                  ***
     **************************************************************************************************
     **************************************************************************************************
     **************************************************************************************************/

    class EigenMatFile {

            /**************************************************************************************************
             **************************************************************************************************
             **                                        PUBLIC SECTION                                        **
             **************************************************************************************************
             **************************************************************************************************/

        public:

            /**************************************************************************************************
             *                                  CONSTRUCTORS AND DESTRUCTOR                                   *
             **************************************************************************************************/


            inline
            EigenMatFile(
                const char* fileName
            );

            // default constructor --> need to .open() it later!
            inline
            EigenMatFile(
            );

            inline
            ~EigenMatFile(
            );


            /**************************************************************************************************
             *                                            READING                                             *
             **************************************************************************************************/

            // read a variable from a mat file;

            template <typename EigenMatrixType>
            inline
            EigenMatrixType
            readVariable(
                const char* variableName
            );


            /**************************************************************************************************
             *                                            WRITING                                             *
             **************************************************************************************************/

            // truncate (delete contents of) mat file
            inline
            void
            truncate(
            );

            // write a variable to a mat file;
            inline
            void
            writeVariable(
                const char* variableName,
                ConstMatrixRef matrix
            );

            inline
            void
            writeVariable(
                const char* variableName,
                ConstSparseMatrixRef matrix
            );

            inline
            void
            writeVariable(
                const char* variableName,
                ConstComplexMatrixRef matrix
            );

            inline
            void
            writeVariable(
                const char* variableName,
                ConstSparseComplexMatrixRef matrix
            );


            inline
            void
            writeVariable(
                const char* variableName,
                const Scalar value
            );

            inline
            void
            writeVariable(
                const char* variableName,
                const ComplexScalar value
            );


            /**************************************************************************************************
             *                                  EXPLICIT OPENING AND CLOSING                                  *
             **************************************************************************************************/

            inline
            void
            open(
                const char* fileName
            );

            inline
            void
            close(
            );


            /**************************************************************************************************
             **************************************************************************************************
             **                                       PRIVATE SECTION                                        **
             **************************************************************************************************
             **************************************************************************************************/

        private:
            char mode_M;
            const char* filename_M;
            SMIO_MatFileHandle_t matFile_M;
            std::unordered_set<std::string> availableVariables_M;

            inline
            void
            reopen(
                const char newMode
            );

            // look up all variables which are there in the mat file.
            // This enables updates of variables in a mat file
            inline
            void
            populateAvailableVariables(
            );

            template <typename EigenMatrixType>
            static inline
            void
            checkDimensions(
                const std::size_t rowsAtRunTime,
                const std::size_t colsAtRunTime
            );

            // actually read the matrix
            // dense and real
            template <typename EigenMatrixType>
            inline
            std::enable_if_t<
                SMIO::isRealMatrix<EigenMatrixType>::value && !SMIO::isSparseMatrix<EigenMatrixType>::value,
                EigenMatrixType
            >
            internalReadVariable(
                const char* variableName
            );

            // dense and complex
            template <typename EigenMatrixType>
            inline
            std::enable_if_t<
                !SMIO::isRealMatrix<EigenMatrixType>::value && !SMIO::isSparseMatrix<EigenMatrixType>::value,
                EigenMatrixType
            >
            internalReadVariable(
                const char* variableName
            );

            // sparse and real
            template <typename EigenMatrixType>
            inline
            std::enable_if_t<
                SMIO::isRealMatrix<EigenMatrixType>::value && SMIO::isSparseMatrix<EigenMatrixType>::value,
                EigenMatrixType
            >
            internalReadVariable(
                const char* variableName
            );

            // sparse and complex
            template <typename EigenMatrixType>
            inline
            std::enable_if_t<
                !SMIO::isRealMatrix<EigenMatrixType>::value && SMIO::isSparseMatrix<EigenMatrixType>::value,
                EigenMatrixType
            >
            internalReadVariable(
                const char* variableName
            );

            // scalar
            template <typename EigenMatrixType>
            inline
            std::enable_if_t<
                std::is_same<EigenMatrixType, SMIO::Scalar>::value || std::is_same<EigenMatrixType, SMIO::ComplexScalar>::value,
                EigenMatrixType
            >
            internalReadVariable(
                const char* variableName
            );

            inline
            bool
            deleteIfAlreadyPresent(
                const std::string& variableName
            );

            inline
            void
            assertReadyForReading(
                const std::string& variableName
            );

            inline
            void
            assertReadyForWriting(
                const std::string& variableName
            );

            [[noreturn]] static inline
            void
            rethrowSMIOError(
                const char* errorMessage
            );

    };


}


#include "EigenSimplematio.cpp"

