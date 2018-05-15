#pragma once
#include "threadLocal.h"
#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include "simplematio.h"
#include "noReturn.h"

#if defined SMIO_DEBUG
#define DEBUG_STATEMENT(command) do {command;} while (false);
#else
#define DEBUG_STATEMENT(command)
#endif

/**************************************************************************************************
 **************************************************************************************************
 **                                      INTERNAL VARIABLES                                      **
 **************************************************************************************************
 **************************************************************************************************/

/**************************************************************************************************
 *                                     SCALAR TYPE SWITCHING                                      *
 **************************************************************************************************/

#ifdef SMIO_USE_FLOATS
SMIO_PRIVATE_API extern
const enum matio_classes SMIO_MatioDenseClass;
SMIO_PRIVATE_API extern
const enum matio_types SMIO_MatioType;

#ifdef CMPLXF
#define MAKE_COMPLEX(x, y) CMPLXF(x, y)
#elif (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) && !defined(__INTEL_COMPILER)
#define MAKE_COMPLEX(x, y) __builtin_complex((float) x, (float) y)
#else
#define MAKE_COMPLEX(x, y) ((float _Complex) ((x) + (y) * _Complex_I))
#endif


#define GET_REAL_PART(x) crealf(x)
#define GET_IMAGINARY_PART(x) cimagf(x)
#define FORMATCODE_FOR_UNDERLYING_DATATYPE "%f"
#define UNDERLYING_DATATYPE_AS_STRING "float"



#else
SMIO_PRIVATE_API extern
const enum matio_classes SMIO_MatioDenseClass;
SMIO_PRIVATE_API extern
const enum matio_types SMIO_MatioType;

#ifdef CMPLXF
#define MAKE_COMPLEX(x, y) CMPLX(x, y)
#elif (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)) && !defined(__INTEL_COMPILER)
#define MAKE_COMPLEX(x, y) __builtin_complex((double) x, (double) y)
#else
#define MAKE_COMPLEX(x, y) ((double _Complex) ((x) + (y) * _Complex_I))
#endif


#define GET_REAL_PART(x) creal(x)
#define GET_IMAGINARY_PART(x) cimag(x)
#define FORMATCODE_FOR_UNDERLYING_DATATYPE "%lf"
#define UNDERLYING_DATATYPE_AS_STRING "double"
#endif

/**************************************************************************************************
 *                                        OTHER VARIABLES                                         *
 **************************************************************************************************/

SMIO_PRIVATE_API extern
const enum matio_classes SMIO_MatioSparseClass;
SMIO_PRIVATE_API extern
THREAD_LOCAL enum matio_compression SMIO_Compression;


#define SMIO_ERROR_MESSAGE_BUFFER_SIZE 200
SMIO_PRIVATE_API
THREAD_LOCAL char errorMessageBuffer[SMIO_ERROR_MESSAGE_BUFFER_SIZE];
SMIO_PRIVATE_API extern
THREAD_LOCAL SMIO_ErrorHandler_t errorHandler;

/**************************************************************************************************
 **************************************************************************************************
 **                                        ERROR HANDLING                                        **
 **************************************************************************************************
 **************************************************************************************************/

SMIO_PRIVATE_API NO_RETURN
void
SMIO_DefaultErrorHandler(
    const char* errorMessage
);

SMIO_PRIVATE_API
void
SMIO_Fail(
    const char* formatString,
    ...
);

/**************************************************************************************************
 **************************************************************************************************
 **                                           READING                                            **
 **************************************************************************************************
 **************************************************************************************************/

SMIO_PRIVATE_API
matvar_t*
SMIO_ReadVariableOrDie(
    mat_t* mat,
    const char* variableName
);

/**************************************************************************************************
 **************************************************************************************************
 **                                           WRITING                                            **
 **************************************************************************************************
 **************************************************************************************************/

SMIO_PRIVATE_API
void
SMIO_FailIfDenseMatrixTooBig(
    const char* variableName,
    const size_t rows,
    const size_t cols,
    const size_t bytesPerEntry
);

SMIO_PRIVATE_API
void
SMIO_FailIfSparseMatrixTooBig(
    const char* variableName,
    const size_t numberOfNonzeroEntries,
    const size_t bytesPerEntry
);

/**************************************************************************************************
 **************************************************************************************************
 **                                        MATRIX CHECKS                                         **
 **************************************************************************************************
 **************************************************************************************************/

SMIO_PRIVATE_API
void
SMIO_AssertIsThatTypeOfMatrix(
    matvar_t* variable,
    const bool shouldBeComplex,
    const bool shouldBeSparse
);

SMIO_PRIVATE_API
void
SMIO_AssertIsScalar(
    const size_t rows,
    const size_t cols,
    const char* name
);

/**************************************************************************************************
 **************************************************************************************************
 **                                    COMPLEX NUMBER SUPPORT                                    **
 **************************************************************************************************
 **************************************************************************************************/

SMIO_PRIVATE_API
SMIO_ComplexScalar_t
SMIO_MakeComplex(
    const SMIO_RealScalar_t realPart,
    const SMIO_RealScalar_t imaginaryPart
);

SMIO_PRIVATE_API
void
SMIO_C99ToSplitComplex(
    const SMIO_ComplexScalar_t* c99ComplexNumbers,
    const size_t numberOfEntries,
    SMIO_RealScalar_t** realParts,
    SMIO_RealScalar_t** imaginaryParts
);

SMIO_PRIVATE_API
void
SMIO_SplitToC99Complex(
    SMIO_ComplexScalar_t** c99ComplexNumbers,
    const size_t numberOfEntries,
    const SMIO_RealScalar_t* realParts,
    const SMIO_RealScalar_t* imaginaryParts
);


