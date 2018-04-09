//
// Created by Stefano Weidmann on 30.03.18.
//
#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>
#include <stefCommonHeaders/dev.hpp>

Matrix computeFluxFunctionFactors(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil) {
    const auto map = [&] (const Real saturationWater) {
        const Real saturationOil = 1.0 - saturationWater;
        const Real saturationOilSquared = saturationOil * saturationOil;
        const Real saturationWaterSquared = saturationWater * saturationWater;
        return  saturationWaterSquared / (porosity * (saturationWaterSquared + saturationOilSquared * dynamicViscosityWater / dynamicViscosityOil));
    };

    return saturationsWater.unaryExpr(map);
}

Matrix computeFluxesX(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesX) {
    ASSERT(fluxFunctionFactors.cols() == darcyVelocitiesX.cols() - 1);
    const int numberOfCols = darcyVelocitiesX.cols();
    for (int colIndex = 1; colIndex < numberOfCols - 1; ++colIndex) {
        darcyVelocitiesX.col(colIndex) *= (darcyVelocitiesX.col(colIndex).array() >= 0).select(
              fluxFunctionFactors.col(colIndex - 1), fluxFunctionFactors.col(colIndex));
    }
    return darcyVelocitiesX;
}

Matrix computeFluxesY(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesY) {
    ASSERT(fluxFunctionFactors.rows() == darcyVelocitiesY.rows() - 1);
    const int numberOfRows = darcyVelocitiesY.rows();
    for (int rowIndex = 1; rowIndex < numberOfRows - 1; ++rowIndex) {
        darcyVelocitiesY.row(rowIndex) *= (darcyVelocitiesY.row(rowIndex).array() >= 0).select(
              fluxFunctionFactors.row(rowIndex), fluxFunctionFactors.row(rowIndex - 1));
    }
    return darcyVelocitiesY;
}


SparseMatrix computePressureResidualsDerivedByPressure(const SparseMatrix& transmissibilities) {
    return transmissibilities;
}

Matrix computePressureResidualsDerivedBySaturationWater() {
    #pragma message "TODO"
    DEV_STUB();
}

Matrix computeSaturationWaterResidualsDerivedBySaturationWater() {
    #pragma message "TODO"
    DEV_STUB();
}

Matrix computeSaturationWaterResidualsDerivedByPressure() {
    #pragma message "TODO"
    DEV_STUB();
}



Matrix advanceStateInTime(Matrix state, ConstMatrixRef derivativeInTime, const Real timestep) {
    state += timestep * derivativeInTime;
}

Matrix computeSaturationDivergence(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, const Real meshWidth) {
    const int numberOfRows = fluxFunctionFactors.rows();
    const int numberOfCols = fluxFunctionFactors.cols();
    ASSERT(fluxesX.cols() == numberOfCols + 1);
    ASSERT(fluxesX.rows() == numberOfRows);
    ASSERT(fluxesY.cols() == numberOfCols);
    ASSERT(fluxesY.rows() == numberOfRows+1);

    Matrix divergences(numberOfRows, numberOfCols);

    for (int colIndex = 0; colIndex < numberOfCols; ++colIndex) {
        divergences.col(colIndex) = (fluxesX.col(colIndex + 1) - fluxesX.col(colIndex)) / meshWidth;
    }

    for (int rowIndex = 0; rowIndex < numberOfRows; ++rowIndex) {
        divergences.row(rowIndex) = (fluxesY.row(rowIndex) - fluxesY.row(rowIndex+1)) / meshWidth;
    }

    return divergences;
}