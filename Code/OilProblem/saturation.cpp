//
// Created by Stefano Weidmann on 30.03.18.
//
#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>
#include <stefCommonHeaders/dev.hpp>
#include "logging.hpp"
#include <stefCommonHeaders/dbg.hpp>

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
        darcyVelocitiesX.col(colIndex).array() *= (darcyVelocitiesX.col(colIndex).array() >= 0).select(
              fluxFunctionFactors.col(colIndex - 1), fluxFunctionFactors.col(colIndex)).array();
    }
    return darcyVelocitiesX;
}

Matrix computeFluxesY(ConstMatrixRef fluxFunctionFactors, Matrix darcyVelocitiesY) {
    ASSERT(fluxFunctionFactors.rows() == darcyVelocitiesY.rows() - 1);
    const int numberOfRows = darcyVelocitiesY.rows();
    for (int rowIndex = 1; rowIndex < numberOfRows - 1; ++rowIndex) {
        darcyVelocitiesY.row(rowIndex).array() *= (darcyVelocitiesY.row(rowIndex).array() >= 0).select(
              fluxFunctionFactors.row(rowIndex), fluxFunctionFactors.row(rowIndex - 1)).array();
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



Matrix computeSaturationDivergences(ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxesX, ConstMatrixRef fluxesY,
                                    const Real meshWidth) {
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

static inline Real computeCflTimestep(const Real maximumAdvectionVelocity, const Real meshWidth) {
    return 0.5 * meshWidth / maximumAdvectionVelocity;
}

void advanceSaturationsInTime(const FixedParameters& params, MatrixRef saturationsWater,
                              ConstMatrixRef pressures,
                              const Eigen::Ref<const Eigen::Matrix<double, -1, -1>>& totalMobilities, Real& time) {

    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);

    LOGGER->debug("pressures {}", pressures);
    const Matrix pressureDerivativesX = computeXDerivative(pressures, params.meshWidth);
    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix pressureDerivativesY = computeYDerivative(pressures, params.meshWidth);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);

    const Matrix fluxesX = computeFluxesX(fluxFunctionFactors, darcyVelocitiesX);
    const Matrix fluxesY = computeFluxesY(fluxFunctionFactors, darcyVelocitiesY);

    const Matrix advectionVelocitiesX = computeXDerivative(fluxFunctionFactors, params.meshWidth).cwiseProduct(darcyVelocitiesX);
    const Matrix advectionVelocitiesY = computeYDerivative(fluxFunctionFactors, params.meshWidth).cwiseProduct(darcyVelocitiesY);

    LOGGER->debug("advection velocities x {}", advectionVelocitiesX);

    const Real maximumAdvectionSpeedX = advectionVelocitiesX.middleCols(1, advectionVelocitiesX.cols()-1).cwiseAbs().maxCoeff();
    const Real maximumAdvectionSpeedY = advectionVelocitiesY.middleRows(1, advectionVelocitiesY.rows()-1).cwiseAbs().maxCoeff();


    const Real timestepX = computeCflTimestep(maximumAdvectionSpeedX, params.meshWidth);
    const Real timestepY = computeCflTimestep(maximumAdvectionSpeedY, params.meshWidth);

    LOGGER->debug("max advection speed x = {}", maximumAdvectionSpeedX);
    LOGGER->debug("max advection speed y = {}", maximumAdvectionSpeedY);

    const Real timestep = std::min(timestepX, timestepY);
    LOGGER->debug("CFL timestep = {}", timestep);

    const Matrix saturationDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);

    const CellIndex drillCell = {saturationsWater.rows()-1, 0};
    const CellIndex wellCell = {0, saturationsWater.cols()-1};

    saturationsWater -= timestep * saturationDivergences;
    wellCell(saturationsWater) -= timestep * params.outflowPerUnitDepthWater(time);
    drillCell(saturationsWater) += timestep * params.inflowPerUnitDepthWater(time);

    time += timestep;

}