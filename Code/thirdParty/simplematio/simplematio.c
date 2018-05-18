#include "simplematioPrivate.h"



/**************************************************************************************************
 **************************************************************************************************
 **                                      INTERNAL VARIABLES                                      **
 **************************************************************************************************
 **************************************************************************************************/

/**************************************************************************************************
 *                                     SCALAR TYPE SWITCHING                                      *
 **************************************************************************************************/

#ifdef SMIO_USE_FLOATS
const enum matio_classes SMIO_MatioDenseClass = MAT_C_SINGLE;
const enum matio_types SMIO_MatioType = MAT_T_SINGLE;

#else
const enum matio_classes SMIO_MatioDenseClass = MAT_C_DOUBLE;
const enum matio_types SMIO_MatioType = MAT_T_DOUBLE;
#endif

/**************************************************************************************************
 *                                        OTHER VARIABLES                                         *
 **************************************************************************************************/

const enum matio_classes SMIO_MatioSparseClass = MAT_C_SPARSE;
THREAD_LOCAL enum matio_compression SMIO_Compression = MAT_COMPRESSION_ZLIB;


THREAD_LOCAL SMIO_ErrorHandler_t errorHandler = SMIO_DefaultErrorHandler;

/**************************************************************************************************
 **************************************************************************************************
 **                                        ERROR HANDLING                                        **
 **************************************************************************************************
 **************************************************************************************************/

void
SMIO_DefaultErrorHandler(
    const char* errorMessage
) {
    fprintf(stderr, "Simple MatIO: ");
    fputs(errorMessage, stderr);
    fputs("", stderr);
    abort();
}


void
SMIO_Fail(
    const char* formatString,
    ...
) {
    va_list argsList;
    va_start(argsList, formatString);

    vsnprintf(errorMessageBuffer, SMIO_ERROR_MESSAGE_BUFFER_SIZE, formatString, argsList);
    // vsnprintf maybe doesn't write \0
    errorMessageBuffer[SMIO_ERROR_MESSAGE_BUFFER_SIZE - 1] = '\0';

    va_end(argsList);

    errorHandler(errorMessageBuffer);
}

SMIO_ErrorHandler_t
SMIO_GetErrorHandler(
) {
    return errorHandler;
}

void
SMIO_SetErrorHandler(
    SMIO_ErrorHandler_t newErrorHandler
) {
    errorHandler = newErrorHandler;
}



/**************************************************************************************************
 **************************************************************************************************
 **                                        MATRIX CHECKS                                         **
 **************************************************************************************************
 **************************************************************************************************/

void
SMIO_AssertIsThatTypeOfMatrix(
    matvar_t* variable,
    const bool shouldBeComplex,
    const bool shouldBeSparse
) {
    const char* name = variable->name;
    const size_t rank = variable->rank;
    if (!shouldBeSparse && variable->class_type != SMIO_MatioDenseClass) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Datatype of variable \"%s\" isn't " UNDERLYING_DATATYPE_AS_STRING ", but this version is configured for " UNDERLYING_DATATYPE_AS_STRING " usage only!",
            name
        );
    } else if (shouldBeSparse && variable->class_type != SMIO_MatioSparseClass) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Variable \"%s\" isn't sparse but it should!",
            name
        );
    } else if (variable->isLogical) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Variable \"%s\" is a logical, but it shouldn't!",
            name
        );
    } else if (variable->isComplex && !shouldBeComplex) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Variable \"%s\" is complex, but should be real!",
            name
        );
    } else if (!variable->isComplex && shouldBeComplex) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Variable \"%s\" is real, but should be complex!",
            name
        );
    } else if (rank != 2) {

        Mat_VarFree(variable);
        SMIO_Fail(
            "Variable \"%s\" has %zu dimensions instead of 2!",
            name,
            rank
        );
    }

    for (size_t dim = 0; dim < 2; ++dim) {
        if (variable->dims[dim] == 0) {
            SMIO_Fail(
                "Variable \"%s\" has zero length in dimension %zu (>= 0)!",
                name,
                dim
            );
        }
    }

    return;
}

void
SMIO_AssertIsScalar(
    const size_t rows,
    const size_t cols,
    const char* name
) {

    const size_t dims[2] = {rows, cols};
    const static char* dimNames[2] = {"rows", "cols"};
    
    for (size_t dim = 0; dim < 2; ++dim) {
        const size_t sizeOfDimension = dims[dim];
        if (sizeOfDimension != 1) {
            SMIO_Fail(
                "Variable \"%s\" has %zu %s but it's supposed to be a scalar!",
                name,
                sizeOfDimension,
                dimNames[dim]
            );
        }
    }
}





/**************************************************************************************************
 **************************************************************************************************
 **                                    COMPLEX NUMBER SUPPORT                                    **
 **************************************************************************************************
 **************************************************************************************************/

SMIO_ComplexScalar_t
SMIO_MakeComplex(
    const SMIO_RealScalar_t realPart,
    const SMIO_RealScalar_t imaginaryPart
) {
    const SMIO_ComplexScalar_t z = MAKE_COMPLEX(realPart, imaginaryPart);
    assert(GET_REAL_PART(z) == realPart && GET_IMAGINARY_PART(z) == imaginaryPart);
    return z;
}

void
SMIO_C99ToSplitComplex(
    const SMIO_ComplexScalar_t* c99ComplexNumbers,
    const size_t numberOfEntries,
    SMIO_RealScalar_t** realParts,
    SMIO_RealScalar_t** imaginaryParts
) {
    assert(realParts && imaginaryParts && c99ComplexNumbers);

    *realParts = malloc(numberOfEntries * sizeof(SMIO_RealScalar_t));
    *imaginaryParts = malloc(numberOfEntries * sizeof(SMIO_RealScalar_t));
    if (*realParts == NULL || *imaginaryParts == NULL) {
        SMIO_Fail(
            "Failed to convert %lu C99 complex numbers to split complex; too little memory available!",
            numberOfEntries
        );
    }

    SMIO_RealScalar_t* x = *realParts;
    SMIO_RealScalar_t* y = *imaginaryParts;
    const SMIO_ComplexScalar_t* z = c99ComplexNumbers;
    for (size_t entry = 0; entry < numberOfEntries; ++entry, ++x, ++y, ++z) {
        *x = GET_REAL_PART(*z);
        *y = GET_IMAGINARY_PART(*z);
    }

    return;
}

void
SMIO_SplitToC99Complex(
    SMIO_ComplexScalar_t** c99ComplexNumbers,
    const size_t numberOfEntries,
    const SMIO_RealScalar_t* realParts,
    const SMIO_RealScalar_t* imaginaryParts
) {
    assert(realParts && imaginaryParts && c99ComplexNumbers);
    *c99ComplexNumbers = malloc(numberOfEntries * sizeof(SMIO_ComplexScalar_t));
    if (*c99ComplexNumbers == NULL) {
        SMIO_Fail(
            "Failed to convert %lu split complex numbers to C99; too little memory available!",
            numberOfEntries
        );
    }

    const SMIO_RealScalar_t* x = realParts;
    const SMIO_RealScalar_t* y = imaginaryParts;
    SMIO_ComplexScalar_t* z = *c99ComplexNumbers;

    for (size_t entry = 0; entry < numberOfEntries; ++entry, ++x, ++y, ++z) {
        *z = SMIO_MakeComplex(*x, *y);
    }

    return;
}

/**************************************************************************************************
 **************************************************************************************************
 **                                         COMPRESSION                                          **
 **************************************************************************************************
 **************************************************************************************************/


void
SMIO_EnableCompression() {
    SMIO_Compression = MAT_COMPRESSION_ZLIB;
    return;
}

void
SMIO_DisableCompression() {
    SMIO_Compression = MAT_COMPRESSION_NONE;
    return;
}


/**************************************************************************************************
 **************************************************************************************************
 **                                     OPENING AND CLOSING                                      **
 **************************************************************************************************
 **************************************************************************************************/

SMIO_MatFileHandle_t
SMIO_Open(
    const char* fileName,
    const char mode
) {
    SMIO_MatFileHandle_t matFile = NULL;
    switch (mode) {
        case 'r': {
            matFile = Mat_Open(fileName, MAT_ACC_RDONLY);
            break;
        }

        case 'w': {
            matFile = Mat_CreateVer(fileName, NULL, MAT_FT_MAT5);
            break;
        }

        case 'a': {
            matFile = Mat_Open(fileName, MAT_ACC_RDWR);
            break;
        }

        default: {
            SMIO_Fail("Invalid access specifier '%c'. Valid are 'r', 'a', 'w'!", mode);
        }
    }

    if (matFile == NULL) {
        SMIO_Fail("Couldn't open file \"%s\" with access specifier '%c'!", fileName, mode);
    }
    return matFile;
}

void
SMIO_Close(
    SMIO_MatFileHandle_t matFile
) {
    assert(matFile);
    Mat_Close(matFile);
    return;
}


/**************************************************************************************************
 **************************************************************************************************
 **                                           READING                                            **
 **************************************************************************************************
 **************************************************************************************************/

matvar_t*
SMIO_ReadVariableOrDie(
    mat_t* mat,
    const char* variableName
) {
    matvar_t* variable = Mat_VarRead(mat, variableName);
    if (variable == NULL) {
        SMIO_Fail(
            "Error reading variable \"%s\" from mat file \"%s\"! Does the variable even exist?",
            variableName,
            Mat_GetFilename(mat)
        );
    }

    return variable;
}

SMIO_DenseRealMatrix_t
SMIO_ReadDenseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
) {
    assert(matFile);

    matvar_t* variable = SMIO_ReadVariableOrDie(matFile, variableName);
    SMIO_AssertIsThatTypeOfMatrix(variable, false /* not complex */, false /* not sparse */);

    SMIO_DenseRealMatrix_t realMatrix = {
        .rows = variable->dims[0],
        .cols = variable->dims[1],
        .data = (SMIO_RealScalar_t*) variable->data
    };

    variable->data = NULL; // Trick to prevent Mat_VarFree from deleting the data
    Mat_VarFree(variable);

    return realMatrix;
}

SMIO_DenseComplexMatrix_t
SMIO_ReadDenseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
) {
    assert(matFile);

    matvar_t* variable = SMIO_ReadVariableOrDie(matFile, variableName);
    SMIO_AssertIsThatTypeOfMatrix(variable, true /* complex */, false /* not sparse */);

    DEBUG_STATEMENT(Mat_VarPrint(variable, 1));

    SMIO_DenseComplexMatrix_t complexMatrix = {
        .rows = variable->dims[0],
        .cols = variable->dims[1],
        .data = NULL
    };

    mat_complex_split_t* splitComplexData = (mat_complex_split_t*) variable->data;
    const SMIO_RealScalar_t* realPart = splitComplexData->Re;
    const SMIO_RealScalar_t* imaginaryPart = splitComplexData->Im;
    const size_t numberOfEntries = complexMatrix.rows * complexMatrix.cols;

    SMIO_ComplexScalar_t* c99ComplexNumbers;

    SMIO_SplitToC99Complex(
        &c99ComplexNumbers,
        numberOfEntries,
        realPart,
        imaginaryPart
    );

    complexMatrix.data = c99ComplexNumbers;

    Mat_VarFree(variable);
    return complexMatrix;
}

SMIO_SparseRealMatrix_t
SMIO_ReadSparseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
) {
    assert(matFile);
    matvar_t* variable = SMIO_ReadVariableOrDie(matFile, variableName);
    if (variable == NULL) {
        exit(5);
    }
    SMIO_AssertIsThatTypeOfMatrix(variable, false /* not complex */, true /* sparse */);


    mat_sparse_t* matioSparse = (mat_sparse_t*) variable->data;

    SMIO_SparseRealMatrix_t matrix = {
        .rows = variable->dims[0],
        .cols = variable->dims[1],
        .data = matioSparse->data,
        .numberOfNonzeroEntries = matioSparse->ndata,
        .rowIndices = matioSparse->ir,
        .numberOfRowIndices = matioSparse->nir,
        .colStarts = matioSparse->jc,
        .numberOfColStarts = matioSparse->njc
    };


    DEBUG_STATEMENT(puts("Read sparse variable:"));
    DEBUG_STATEMENT(printf("ndata = %d, nir = %d, njc = %d\n", matioSparse->ndata, matioSparse->nir, matioSparse->njc));
    DEBUG_STATEMENT(Mat_VarPrint(variable, 1));

    // Prevent MATIO from deleting the arrays
    variable->data = NULL;


    Mat_VarFree(variable);
    return matrix;
}

SMIO_SparseComplexMatrix_t
SMIO_ReadSparseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
) {
    assert(matFile);
    matvar_t* variable = SMIO_ReadVariableOrDie(matFile, variableName);
    SMIO_AssertIsThatTypeOfMatrix(variable, true /* complex */, true /* sparse */);


    mat_sparse_t* matioSparse = (mat_sparse_t*) variable->data;

    SMIO_SparseComplexMatrix_t matrix = {
        .rows = variable->dims[0],
        .cols = variable->dims[1],
        .data = NULL,
        .numberOfNonzeroEntries = matioSparse->ndata,
        .rowIndices = matioSparse->ir,
        .numberOfRowIndices = matioSparse->nir,
        .colStarts = matioSparse->jc,
        .numberOfColStarts = matioSparse->njc
    };

    // Convert them to C99 complex
    mat_complex_split_t* splitComplexData = (mat_complex_split_t*) matioSparse->data;
    const SMIO_RealScalar_t* realPart = splitComplexData->Re;
    const SMIO_RealScalar_t* imaginaryPart = splitComplexData->Im;

    SMIO_SplitToC99Complex(
        &(matrix.data),
        matrix.numberOfNonzeroEntries,
        realPart,
        imaginaryPart
    );

    // Prevent MATIO from deleting the arrays
    matioSparse->ir = NULL;
    matioSparse->jc = NULL;


    Mat_VarFree(variable);
    return matrix;
}



SMIO_RealScalar_t
SMIO_ReadRealScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
) {
    SMIO_DenseRealMatrix_t scalarAsMatrix = SMIO_ReadDenseRealMatrix(matFile, variableName);
    SMIO_AssertIsScalar(scalarAsMatrix.rows, scalarAsMatrix.cols, variableName);

    return *scalarAsMatrix.data;
}


SMIO_ComplexScalar_t
SMIO_ReadComplexScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName
) {
    SMIO_DenseComplexMatrix_t scalarAsMatrix = SMIO_ReadDenseComplexMatrix(matFile, variableName);
    SMIO_AssertIsScalar(scalarAsMatrix.rows, scalarAsMatrix.cols, variableName);

    return *scalarAsMatrix.data;
}




/**************************************************************************************************
 **************************************************************************************************
 **                                           WRITING                                            **
 **************************************************************************************************
 **************************************************************************************************/

void
SMIO_FailIfDenseMatrixTooBig(
    const char* variableName,
    const size_t rows,
    const size_t cols,
    const size_t bytesPerEntry
) {
    const size_t bytesToWrite = rows * cols * bytesPerEntry;
    const size_t maxBytesToWrite = 1 << 31;  // 2^31
    if (bytesToWrite > maxBytesToWrite) {
        SMIO_Fail("Couldn't write dense matrix %s: can only write %lu bytes at most, but should write %lu!",
                  variableName, maxBytesToWrite, bytesToWrite);
    }
}

void
SMIO_FailIfSparseMatrixTooBig(
    const char* variableName,
    const size_t numberOfNonzeroEntries,
    const size_t bytesPerEntry
) {
    const size_t bytesToWrite = numberOfNonzeroEntries* bytesPerEntry;
    const size_t maxBytesToWrite = 1 << 31;  // 2^31
    if (bytesToWrite > maxBytesToWrite) {
        SMIO_Fail("Couldn't write sparse matrix %s: can only write %lu bytes at most, but should write %lu!",
                  variableName, maxBytesToWrite, bytesToWrite);
    }
}


void
SMIO_WriteDenseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_DenseRealMatrix_t* matrix
) {
    assert(matFile && matrix->data);

    SMIO_FailIfDenseMatrixTooBig(variableName, matrix->rows, matrix->cols, sizeof(SMIO_RealScalar_t));

    size_t dimensions[2] = {matrix->rows, matrix->cols};
    int options = MAT_F_DONT_COPY_DATA;

    matvar_t* variable = Mat_VarCreate(
                             variableName,
                             SMIO_MatioDenseClass,
                             SMIO_MatioType,
                             2, /* a matrix has 2 dimensions */
                             dimensions,
                             matrix->data,
                             options
                         );

    DEBUG_STATEMENT(printf("Created variable \"%s\":\n", variableName));
    DEBUG_STATEMENT(Mat_VarPrint(variable, 1));

    if (variable == NULL) {
        SMIO_Fail(
            "Couldn't create variable \"%s\", a %lux%lu dense real " UNDERLYING_DATATYPE_AS_STRING " matrix!",
            variableName,
            matrix->rows,
            matrix->cols
        );
    }

    if (Mat_VarWrite(matFile, variable, SMIO_Compression) != 0) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Couldn't write variable \"%s\", a %lux%lu dense real " UNDERLYING_DATATYPE_AS_STRING " matrix, to file \"%s\"!",
            variableName,
            matrix->rows,
            matrix->cols,
            Mat_GetFilename(matFile)
        );
    }


    Mat_VarFree(variable);

    return;
}

void
SMIO_WriteDenseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_DenseComplexMatrix_t* matrix
) {
    assert(matrix->data);

    SMIO_FailIfDenseMatrixTooBig(variableName, matrix->rows, matrix->cols, sizeof(SMIO_ComplexScalar_t));

    size_t dimensions[2] = {matrix->rows, matrix->cols};
    const int options = MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX;

    const size_t numberOfEntries = matrix->rows * matrix->cols;

    SMIO_RealScalar_t* imaginaryParts;
    SMIO_RealScalar_t* realParts;

    SMIO_C99ToSplitComplex(
        matrix->data,
        numberOfEntries,
        &realParts,
        &imaginaryParts
    );

    mat_complex_split_t splitComplexData = {
        .Re = (void*) realParts,
        .Im = (void*) imaginaryParts
    };

    matvar_t* variable = Mat_VarCreate(
                             variableName,
                             SMIO_MatioDenseClass,
                             SMIO_MatioType,
                             2,
                             dimensions,
                             (void*) &splitComplexData,
                             options
                         );

    DEBUG_STATEMENT(printf("Created variable \"%s\":\n", variableName));
    DEBUG_STATEMENT(Mat_VarPrint(variable, 1));

    if (variable == NULL) {
        SMIO_Fail(
            "Couldn't create variable \"%s\", a %lux%lu dense complex " UNDERLYING_DATATYPE_AS_STRING " matrix!",
            variableName,
            matrix->rows,
            matrix->cols
        );
    }

    if (Mat_VarWrite(matFile, variable, SMIO_Compression) != 0) {
        Mat_VarFree(variable);
        free(imaginaryParts);
        free(realParts);
        SMIO_Fail(
            "Couldn't write variable \"%s\", a %lux%lu dense complex " UNDERLYING_DATATYPE_AS_STRING " matrix, to file \"%s\"!",
            variableName,
            matrix->rows,
            matrix->cols,
            Mat_GetFilename(matFile)
        );
    }

    Mat_VarFree(variable);
    free(imaginaryParts);
    free(realParts);

    return;
}

void
SMIO_WriteSparseRealMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_SparseRealMatrix_t* matrix
) {

    assert(matFile && matrix);
    SMIO_FailIfSparseMatrixTooBig(variableName, matrix->numberOfNonzeroEntries, sizeof(SMIO_RealScalar_t));

    size_t dimensions[2] = {matrix->rows, matrix->cols};
    int options = MAT_F_DONT_COPY_DATA;

    mat_sparse_t matioSparse = {
        .nzmax = matrix->numberOfNonzeroEntries,
        .ir = matrix->rowIndices,
        .nir = matrix->numberOfRowIndices,
        .jc = matrix->colStarts,
        .njc = matrix->numberOfColStarts,
        .ndata = matrix->numberOfNonzeroEntries,
        .data = matrix->data
    };

    matvar_t* variable = Mat_VarCreate(
                             variableName,
                             SMIO_MatioSparseClass,
                             SMIO_MatioType,
                             2,
                             dimensions,
                             (void*) &matioSparse,
                             options
                         );


    if (variable == NULL) {
        SMIO_Fail(
            "Couldn't create variable \"%s\", a %lux%lu sparse real " UNDERLYING_DATATYPE_AS_STRING " matrix!",
            variableName,
            matrix->rows,
            matrix->cols
        );
    }

    DEBUG_STATEMENT(printf("Created sparse variable \"%s\":\n", variableName));
    DEBUG_STATEMENT(printf("ndata = %d, nir = %d, njc = %d\n", matioSparse.ndata, matioSparse.nir, matioSparse.njc));
    DEBUG_STATEMENT(Mat_VarPrint(variable, 1));

    if (Mat_VarWrite(matFile, variable, SMIO_Compression) != 0) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Couldn't write variable \"%s\", a %lux%lu sparse real " UNDERLYING_DATATYPE_AS_STRING " matrix, to file \"%s\"!",
            variableName,
            matrix->rows,
            matrix->cols,
            Mat_GetFilename(matFile)
        );
    }

    // Prevent MATIO from freeing the data
    variable->data = NULL;
    Mat_VarFree(variable);



    return;
}

void
SMIO_WriteSparseComplexMatrix(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    const SMIO_SparseComplexMatrix_t* matrix
) {

    assert(matFile && matrix && matrix->data);
    SMIO_FailIfSparseMatrixTooBig(variableName, matrix->numberOfNonzeroEntries, sizeof(SMIO_ComplexScalar_t));

    size_t dimensions[2] = {matrix->rows, matrix->cols};
    int options = MAT_F_DONT_COPY_DATA | MAT_F_COMPLEX;

    SMIO_RealScalar_t* imaginaryParts;
    SMIO_RealScalar_t* realParts;

    SMIO_C99ToSplitComplex(
        matrix->data,
        matrix->numberOfNonzeroEntries,
        &realParts,
        &imaginaryParts
    );

    mat_complex_split_t splitComplexData = {
        .Re = (void*) realParts,
        .Im = (void*) imaginaryParts
    };

    mat_sparse_t matioSparse = {
        .nzmax = matrix->numberOfNonzeroEntries,
        .ir = matrix->rowIndices,
        .nir = matrix->numberOfRowIndices,
        .jc = matrix->colStarts,
        .njc = matrix->numberOfColStarts,
        .ndata = matrix->numberOfNonzeroEntries,
        .data = &splitComplexData
    };

    matvar_t* variable = Mat_VarCreate(
                             variableName,
                             SMIO_MatioSparseClass,
                             SMIO_MatioType,
                             2,
                             dimensions,
                             (void*) &matioSparse,
                             options
                         );


    if (variable == NULL) {
        SMIO_Fail(
            "Couldn't create variable \"%s\", a %lux%lu sparse complex " UNDERLYING_DATATYPE_AS_STRING " matrix!",
            variableName,
            matrix->rows,
            matrix->cols
        );
    }

    DEBUG_STATEMENT(printf("Created sparse complex variable \"%s\":\n", variableName));
    DEBUG_STATEMENT(printf("ndata = %d, nir = %d, njc = %d\n", matioSparse.ndata, matioSparse.nir, matioSparse.njc));
    DEBUG_STATEMENT(Mat_VarPrint(variable, 1));

    if (Mat_VarWrite(matFile, variable, SMIO_Compression) != 0) {
        Mat_VarFree(variable);
        SMIO_Fail(
            "Couldn't write variable \"%s\", a %lux%lu sparse complex " UNDERLYING_DATATYPE_AS_STRING " matrix, to file \"%s\"!",
            variableName,
            matrix->rows,
            matrix->cols,
            Mat_GetFilename(matFile)
        );
    }

    // Prevent MATIO from freeing the data
    variable->data = NULL;
    Mat_VarFree(variable);

    free(imaginaryParts);
    free(realParts);

    return;
}


void
SMIO_WriteRealScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    SMIO_RealScalar_t scalar
) {
    SMIO_DenseRealMatrix_t scalarAsMatrix = SMIO_MakeDenseRealMatrix(1, 1, &scalar);
    return SMIO_WriteDenseRealMatrix(matFile, variableName, &scalarAsMatrix);
}

void
SMIO_WriteComplexScalar(
    SMIO_MatFileHandle_t matFile,
    const char* variableName,
    SMIO_ComplexScalar_t scalar
) {
    SMIO_DenseComplexMatrix_t scalarAsMatrix = SMIO_MakeDenseComplexMatrix(1, 1, &scalar);
    return SMIO_WriteDenseComplexMatrix(matFile, variableName, &scalarAsMatrix);
}

/**************************************************************************************************
 *                                        FOR CONVENIENCE                                         *
 **************************************************************************************************/

SMIO_DenseRealMatrix_t
SMIO_MakeDenseRealMatrix(
    size_t rows,
    size_t cols,
    SMIO_RealScalar_t* data
) {
    SMIO_DenseRealMatrix_t matrix = {
        .rows = rows,
        .cols = cols,
        .data = data
    };

    return matrix;
}

SMIO_DenseComplexMatrix_t
SMIO_MakeDenseComplexMatrix(
    size_t rows,
    size_t cols,
    SMIO_ComplexScalar_t* data
) {
    SMIO_DenseComplexMatrix_t matrix = {
        .rows = rows,
        .cols = cols,
        .data = data
    };

    return matrix;
}


