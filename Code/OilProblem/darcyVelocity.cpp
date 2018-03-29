//
// Created by Stefano Weidmann on 29.03.18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>

Matrix computeXDerivative(ConstMatrixRef field, const Real meshWidth) {
    Matrix xDerivative(field.rows(), field.cols() + 1);
    const int numberOfCols = int(xDerivative.cols());
    const int lastColIndex = numberOfCols-1;
    xDerivative.col(0).setZero();

    for (int colIndex = 1; colIndex < lastColIndex; ++colIndex) {
        xDerivative.col(colIndex) = (field.col(colIndex) - field.col(colIndex-1)) / meshWidth;
    }
    xDerivative.col(lastColIndex).setZero();

    return xDerivative;
}

Matrix computeYDerivative(ConstMatrixRef field, const Real meshWidth) {
    Matrix yDerivative(field.rows()+1, field.cols());
    const int numberOfRows = int(yDerivative.rows());
    const int lastRowIndex = numberOfRows-1;
    yDerivative.row(0).setZero();

    for (int rowIndex = 1; rowIndex < lastRowIndex; ++rowIndex) {
        yDerivative.row(rowIndex) = (field.row(rowIndex-1) - field.row(rowIndex)) / meshWidth;
    }
    yDerivative.row(lastRowIndex).setZero();

    return yDerivative;
}

Real getDerivativeAtCellBorder(CellIndex cell,
                               ConstMatrixRef xDerivative, ConstMatrixRef yDerivative,
                               const CellIndex::Direction whichBorder) {

    switch (whichBorder) {
        case CellIndex::Direction::NORTH: {
            return cell(yDerivative);
        }
        case CellIndex::Direction::SOUTH: {
            ++cell.i;
            return cell(yDerivative);
        }
        case CellIndex::Direction::WEST: {
            return cell(xDerivative);
        }
        case CellIndex::Direction::EAST: {
            ++cell.j;
            return cell(xDerivative);
        }
    }
}
