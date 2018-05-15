/*
 ad88888ba  88 88b           d88 88888888ba  88          88888888888
d8"     "8b 88 888b         d888 88      "8b 88          88
Y8,         88 88`8b       d8'88 88      ,8P 88          88
`Y8aaaaa,   88 88 `8b     d8' 88 88aaaaaa8P' 88          88aaaaa
  `"""""8b, 88 88  `8b   d8'  88 88""""""'   88          88"""""
        `8b 88 88   `8b d8'   88 88          88          88
Y8a     a8P 88 88    `888'    88 88          88          88
 "Y88888P"  88 88     `8'     88 88          88888888888 88888888888



88b           d88        db   888888888888 88   ,ad8888ba,
888b         d888       d88b       88      88  d8"'    `"8b
88`8b       d8'88      d8'`8b      88      88 d8'        `8b
88 `8b     d8' 88     d8'  `8b     88      88 88          88
88  `8b   d8'  88    d8YaaaaY8b    88      88 88          88
88   `8b d8'   88   d8""""""""8b   88      88 Y8,        ,8P
88    `888'    88  d8'        `8b  88      88  Y8a.    .a8P
88     `8'     88 d8'          `8b 88      88   `"Y8888Y"'
*/
#pragma once
#ifdef __cplusplus
extern "C" {
#endif


#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>
#include "matio.h"

/**************************************************************************************************
 **************************************************************************************************
 **                                      VISIBILITY CONTROL                                      **
 **************************************************************************************************
 **************************************************************************************************/

// A function or variable declared with SMIO_PUBLIC_API
// is meant for use by everybody. On the other hand, if
// it's declared with SMIO_PRIVATE_API, it's only used internally by this library.
#if defined WIN32 || defined __CYGWIN__
#if defined BUILDING_DLL
#define SMIO_PUBLIC_API __declspec(dllexport)
#else
#define SMIO_PUBLIC_API __declspec(dllimport)
#endif
#define SMIO_PRIVATE_API
#else
#define SMIO_PUBLIC_API __attribute__((visibility("default")))
#define SMIO_PRIVATE_API __attribute__((visibility("hidden")))
#endif


/**************************************************************************************************
 **************************************************************************************************
 **                                           TYPEDEFS                                           **
 **************************************************************************************************
 **************************************************************************************************/

/**************************************************************************************************
 *                                          SCALAR TYPES                                          *
 **************************************************************************************************/

// The underlying Scalar types to use
// defaults to double
// To use floats, define the macro SMIO_USE_FLOATS and link with the library "simplematiof"
// instead of "simplematio"
#if defined SMIO_USE_FLOATS
typedef float SMIO_RealScalar_t;
typedef float _Complex SMIO_ComplexScalar_t;
#else
typedef double SMIO_RealScalar_t;
typedef double _Complex SMIO_ComplexScalar_t;
#endif

/**************************************************************************************************
 *                                       DENSE MATRIX TYPES                                       *
 **************************************************************************************************/

typedef struct {
    size_t rows;
    size_t cols;
    SMIO_RealScalar_t* data; // Stored in column major format
} SMIO_DenseRealMatrix_t;

typedef struct {
    size_t rows;
    size_t cols;
    SMIO_ComplexScalar_t* data; // Stored in column major format
} SMIO_DenseComplexMatrix_t;


/**************************************************************************************************
 *                                      SPARSE MATRIX TYPES                                       *
 **************************************************************************************************/

// A sparse matrix in compressed column storage
// See the article at http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000

typedef struct {
    size_t rows;
    size_t cols;
    SMIO_RealScalar_t* data; // 'val' in the article referred

    size_t numberOfNonzeroEntries;
    int* rowIndices; // row_ind in the article referred
    int numberOfRowIndices;
    int* colStarts; // col_ptr in the article referred
    int numberOfColStarts;
} SMIO_SparseRealMatrix_t;

typedef struct {
    size_t rows;
    size_t cols;
    SMIO_ComplexScalar_t* data;

    size_t numberOfNonzeroEntries;
    int* rowIndices;
    int numberOfRowIndices;
    int* colStarts;
    int numberOfColStarts;
} SMIO_SparseComplexMatrix_t;


/**************************************************************************************************
 *                                          OTHER TYPES                                           *
 **************************************************************************************************/

// The type of an error handler used by this library
typedef void (*SMIO_ErrorHandler_t)(const char* errorMessage);


// The mat file type used by <matio.h>
// Uninteresting, just use it as an opaque handle
typedef mat_t* SMIO_MatFileHandle_t;



/**************************************************************************************************
 **************************************************************************************************
 **                                     OPENING AND CLOSING                                      **
 **************************************************************************************************
 **************************************************************************************************/

// Open a matfile for:
// reading: mode 'r'
// reading and writing: mode 'w' (truncates file)
// reading and appending: mode 'a' (doesn't truncate file)
SMIO_PUBLIC_API
SMIO_MatFileHandle_t
SMIO_Open(
    const char* fileName,
    const char mode
);

SMIO_PUBLIC_API
void
SMIO_Close(
    SMIO_MatFileHandle_t matFile
);


/**************************************************************************************************
 **************************************************************************************************
 **                                           READING                                            **
 **************************************************************************************************
 **************************************************************************************************/

// read various types of matrices from an opened MAT file.
// SMIO_MatFileHandle_t is what SMIO_Open in the section 'OPENING AND CLOSING' returns

SMIO_PUBLIC_API
SMIO_DenseRealMatrix_t
SMIO_ReadDenseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
);

SMIO_PUBLIC_API
SMIO_DenseComplexMatrix_t
SMIO_ReadDenseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
);

SMIO_PUBLIC_API
SMIO_SparseRealMatrix_t
SMIO_ReadSparseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
);

SMIO_PUBLIC_API
SMIO_SparseComplexMatrix_t
SMIO_ReadSparseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
);



// convenience functions for reading scalars
SMIO_PUBLIC_API
SMIO_RealScalar_t
SMIO_ReadRealScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
);

SMIO_PUBLIC_API
SMIO_ComplexScalar_t
SMIO_ReadComplexScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
);




/**************************************************************************************************
 **************************************************************************************************
 **                                           WRITING                                            **
 **************************************************************************************************
 **************************************************************************************************/

// Write a dense real matrix with given rows and columns
// Input must be in colmajor order
SMIO_PUBLIC_API
void
SMIO_WriteDenseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_DenseRealMatrix_t* matrix
);

SMIO_PUBLIC_API
void
SMIO_WriteDenseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_DenseComplexMatrix_t* matrix
);

SMIO_PUBLIC_API
void
SMIO_WriteSparseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_SparseRealMatrix_t* matrix
);

SMIO_PUBLIC_API
void
SMIO_WriteSparseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_SparseComplexMatrix_t* matrix
);

// Convenience functions for writing scalars
SMIO_PUBLIC_API
void
SMIO_WriteRealScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    SMIO_RealScalar_t scalar
);

SMIO_PUBLIC_API
void
SMIO_WriteComplexScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    SMIO_ComplexScalar_t scalar
);

/**************************************************************************************************
 *                                        FOR CONVENIENCE                                         *
 **************************************************************************************************/

SMIO_PUBLIC_API
SMIO_DenseRealMatrix_t
SMIO_MakeDenseRealMatrix(
    size_t rows,
    size_t cols,
    SMIO_RealScalar_t* data
);

SMIO_PUBLIC_API
SMIO_DenseComplexMatrix_t
SMIO_MakeDenseComplexMatrix(
    size_t rows,
    size_t cols,
    SMIO_ComplexScalar_t* data
);

/**************************************************************************************************
 **************************************************************************************************
 **                                         COMPRESSION                                          **
 **************************************************************************************************
 **************************************************************************************************/

// disable or enable zlib compression
// the default is enabled compression

SMIO_PUBLIC_API
void
SMIO_EnableCompression();

SMIO_PUBLIC_API
void
SMIO_DisableCompression();


/**************************************************************************************************
 **************************************************************************************************
 **                                        ERROR HANDLING                                        **
 **************************************************************************************************
 **************************************************************************************************/

// sets a new error handler which this library calls if there is an error.
// It needs to have the signature
// void foobar(char* errorMessage)
SMIO_PUBLIC_API
void
SMIO_SetErrorHandler(
    SMIO_ErrorHandler_t newErrorHandler
);

// gets the current error handler
SMIO_PUBLIC_API
SMIO_ErrorHandler_t
SMIO_GetErrorHandler(
);

#ifdef __cplusplus
}
#endif

#ifndef BUILDING_SMIO
#undef SMIO_PUBLIC_API
#undef SMIO_PRIVATE_API
#endif

