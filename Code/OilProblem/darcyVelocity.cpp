//
// Created by Stefano Weidmann on 29.03.18.
//

#include "darcyVelocity.hpp"
#include "pressure.hpp"
#include <stefCommonHeaders/assert.h>

Matrix computeXDerivative(ConstMatrixRef field, const Real meshWidth) {
    Matrix xDerivative(field.rows(), field.cols() + 1);
    const int numberOfCols = xDerivative.cols();
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
    const int numberOfRows = yDerivative.rows();
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

CellIndex borderIndexToCenterIndex(CellIndex borderIndex, const CellIndex::Direction whichBorder) {
    switch (whichBorder) {
        case CellIndex::Direction::NORTH: {
            return borderIndex;
        }
        case CellIndex::Direction::SOUTH: {
            --borderIndex.i;
            return borderIndex;
        }
        case CellIndex::Direction::WEST: {
            return borderIndex;
        }
        case CellIndex::Direction::EAST: {
            --borderIndex.j;
            return borderIndex;
        }
    }
}

CellIndex centerIndexToBorderIndex(CellIndex centerIndex, const CellIndex::Direction whichBorder) {
    switch (whichBorder) {
        case CellIndex::Direction::NORTH: {
            return centerIndex;
        }
        case CellIndex::Direction::SOUTH: {
            ++centerIndex.i;
            return centerIndex;
        }
        case CellIndex::Direction::WEST: {
            return centerIndex;
        }
        case CellIndex::Direction::EAST: {
            ++centerIndex.j;
            return centerIndex;
        }
    }
}


Matrix computeTotalDarcyVelocitiesX(ConstMatrixRef totalMobilities, Matrix pressureDerivativesX) {
    const int numberOfColsOfDerivatives = pressureDerivativesX.cols();
    const int numberOfRowsOfDerivatives = pressureDerivativesX.rows();
    const int numberOfRowsOfPressures = numberOfRowsOfDerivatives;
    const int numberOfColsOfPressures = numberOfColsOfDerivatives - 1;

    CellIndex borderIndex;

    // skip first and last column which are on the outside of the domain and are zero anyways.
    for (borderIndex.j = 1; borderIndex.j < numberOfColsOfPressures; ++borderIndex.j) {
        for (borderIndex.i = 0; borderIndex.i < numberOfRowsOfPressures; ++borderIndex.i) {
            const CellIndex pressureCellWestOfThisBorder = borderIndexToCenterIndex(borderIndex, CellIndex::Direction::EAST); // This border is to the east of the cell
            const CellIndex pressureCellEastOfThisBorder = borderIndexToCenterIndex(borderIndex, CellIndex::Direction::WEST);

            const Real transmissibility = computeTransmissibility(totalMobilities, pressureCellEastOfThisBorder, pressureCellWestOfThisBorder);

            borderIndex(pressureDerivativesX) *= -transmissibility;
        }
    }

    return pressureDerivativesX;
}

Matrix computeTotalDarcyVelocitiesY(ConstMatrixRef totalMobilities, Matrix pressureDerivativesY) {
    const int numberOfColsOfDerivatives = pressureDerivativesY.cols();
    const int numberOfRowsOfDerivatives = pressureDerivativesY.rows();
    const int numberOfRowsOfPressures = numberOfRowsOfDerivatives - 1; // the zero derivatives at the top and the bottom are added
    const int numberOfColsOfPressures = numberOfColsOfDerivatives;

    CellIndex borderIndex;

    for (borderIndex.j = 0; borderIndex.j < numberOfColsOfPressures; ++borderIndex.j) {
        // skip first and last row which are on the outside of the domain and are zero anyways.
        for (borderIndex.i = 1; borderIndex.i < numberOfRowsOfPressures; ++borderIndex.i) {
            const CellIndex pressureCellNorthOfThisBorder = borderIndexToCenterIndex(borderIndex, CellIndex::Direction::SOUTH); // This border is to the south of the cell
            const CellIndex pressureCellSouthOfThisBorder = borderIndexToCenterIndex(borderIndex, CellIndex::Direction::NORTH);

            const Real transmissibility = computeTransmissibility(totalMobilities, pressureCellNorthOfThisBorder, pressureCellSouthOfThisBorder);

            borderIndex(pressureDerivativesY) *= -transmissibility;
        }
    }

    return pressureDerivativesY;
}