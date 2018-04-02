//
// Created by Stefano Weidmann on 30.03.18.
//
#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>

Matrix computeFluxFunctionFactors(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil) {
    const auto map = [&] (const Real saturationWater) {
        const Real saturationOil = 1.0 - saturationWater;
        const Real saturationOilSquared = saturationOil * saturationOil;
        const Real saturationWaterSquared = saturationWater * saturationWater;
        return  saturationWaterSquared / (porosity * (saturationWaterSquared + saturationOilSquared * dynamicViscosityWater / dynamicViscosityOil));
    };

    return saturationsWater.unaryExpr(map);
}

//Matrix approximateFluxFunctionFactorsAtBordersX(ConstMatrixRef fluxFunctionFactors) {
//    Matrix fluxFunctionFactorsAtBordersX(fluxFunctionFactors.rows(), fluxFunctionFactors.cols()+1);
//
//    fluxFunctionFactorsAtBordersX.col(0).setZero();
//
//    const int numberOfCols = int(fluxFunctionFactorsAtBordersX.cols());
//    for (int colIndex = 1; colIndex < numberOfCols-1; ++colIndex) {
//        fluxFunctionFactorsAtBordersX.col(colIndex) = 0.5 * (fluxFunctionFactorsAtBordersX.col(colIndex) - fluxFunctionFactorsAtBordersX.col(colIndex-1));
//    }
//
//    fluxFunctionFactorsAtBordersX.col(numberOfCols-1).setZero();
//    return fluxFunctionFactorsAtBordersX;
//}
//
//Matrix approximateFluxFunctionFactorsAtBordersY(ConstMatrixRef fluxFunctionFactors) {
//    Matrix fluxFunctionFactorsAtBordersY(fluxFunctionFactors.rows(), fluxFunctionFactors.cols()+1);
//
//    fluxFunctionFactorsAtBordersY.row(0).setZero();
//
//    const int numberOfRows = int(fluxFunctionFactorsAtBordersY.rows());
//    for (int rowIndex = 1; rowIndex < numberOfRows-1; ++rowIndex) {
//        fluxFunctionFactorsAtBordersY.row(rowIndex) = 0.5 * (fluxFunctionFactorsAtBordersY.row(rowIndex) - fluxFunctionFactorsAtBordersY.row(rowIndex-1));
//    }
//
//    fluxFunctionFactorsAtBordersY.col(numberOfRows-1).setZero();
//    return fluxFunctionFactorsAtBordersY;
//}

Matrix computeDivergence(ConstMatrixRef xDerivative, ConstMatrixRef yDerivative) {
    Matrix divergence(xDerivative.rows(), xDerivative.cols()-1);
    ASSERT(divergence.rows() == yDerivative.rows() - 1);
    ASSERT(divergence.cols() == yDerivative.cols());

    const int numberOfColumnsInDivergence = divergence.cols();
    const int numberOfRowsInDivergence = divergence.rows();

    for (int divCol = 0; divCol < numberOfColumnsInDivergence; ++divCol) {
        divergence.col(divCol) =
    }

    for (int divRow = 0; divRow < numberOfRowsInDivergence; ++divRow) {
        divergence.row(divRow) += 0.5 * (yDerivative.row(divRow) + yDerivative.row(divRow + 1));
    }

    return divergence;
}

