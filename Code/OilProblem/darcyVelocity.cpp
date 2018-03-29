//
// Created by Stefano Weidmann on 29.03.18.
//

#include "oilProblem.hpp"
#include "darcyVelocityInternal.hpp"
// Gradient zero at boundaries
Matrix computeGradientComponent(ConstMatrixRef field, const Real meshWidth,
                                const DerivativeDirection direction) {


    Matrix gradientComponent;
    switch (direction) {
        case DerivativeDirection::X: {
            gradientComponent.resize(field.rows(), field.cols()+1);
            const int numberOfCols = int(gradientComponent.cols());
            const int lastColIndex = numberOfCols-1;
            gradientComponent.col(0).setZero();

            for (int colIndex = 1; colIndex < lastColIndex; ++colIndex) {
                gradientComponent.col(colIndex) = (field.col(colIndex) - field.col(colIndex-1)) / meshWidth;
            }
            gradientComponent.col(lastColIndex).setZero();

            break;
        }

        case DerivativeDirection::Y: {
            gradientComponent.resize(field.rows()+1, field.cols());
            const int numberOfRows = int(gradientComponent.rows());
            const int lastRowIndex = numberOfRows-1;
            gradientComponent.row(0).setZero();

            for (int rowIndex = 1; rowIndex < lastRowIndex; ++rowIndex) {
                gradientComponent.row(rowIndex) = (field.row(rowIndex) - field.row(rowIndex+1)) / meshWidth;
            }
            gradientComponent.row(lastRowIndex).setZero();
            break;
        }
    }

    return gradientComponent;

}

