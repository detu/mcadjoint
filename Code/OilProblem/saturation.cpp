//
// Created by Stefano Weidmann on 30.03.18.
//
#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>
#include <stefCommonHeaders/dev.hpp>
#include "logging.hpp"
#include <stefCommonHeaders/dbg.hpp>
#include "specialCells.hpp"

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


Real clamp(const Real x, const Real minVal, const Real maxVal) {
    if (x > maxVal) {
        return maxVal;
    }

    if (x < minVal) {
        return minVal;
    }

    return x;
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
        divergences.row(rowIndex) += (fluxesY.row(rowIndex) - fluxesY.row(rowIndex+1)) / meshWidth;
    }

    return divergences;
}

static inline Real computeCflTimestep(const Real maximumAdvectionVelocity, const Real meshWidth) {
    return 0.9 * meshWidth / maximumAdvectionVelocity;
}

void advanceSaturationsInTime(const FixedParameters& params, MatrixRef saturationsWater,
                              ConstMatrixRef pressures,
                              ConstMatrixRef totalMobilities, Real& time) {

    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);

    Matrix pressureDerivativesX = computeXDerivative(pressures, params.meshWidth);
    Matrix pressureDerivativesY = computeYDerivative(pressures, params.meshWidth);
    //adaptPressureGradientsAtWell(params.inflowPerUnitDepthWater(time), totalMobilities, pressures, pressureDerivativesX, pressureDerivativesY, params.meshWidth);

    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);

    const Matrix fluxesX = computeFluxesX(fluxFunctionFactors, darcyVelocitiesX);
    const Matrix fluxesY = computeFluxesY(fluxFunctionFactors, darcyVelocitiesY);

    const Matrix advectionVelocitiesX = computeXDerivative(fluxFunctionFactors, params.meshWidth).cwiseProduct(darcyVelocitiesX);
    const Matrix advectionVelocitiesY = computeYDerivative(fluxFunctionFactors, params.meshWidth).cwiseProduct(darcyVelocitiesY);

    const Real maximumAdvectionSpeedX = advectionVelocitiesX.cwiseAbs().maxCoeff();
    const Real maximumAdvectionSpeedY = advectionVelocitiesY.cwiseAbs().maxCoeff();


    const Real timestepX = computeCflTimestep(maximumAdvectionSpeedX, params.meshWidth);
    const Real timestepY = computeCflTimestep(maximumAdvectionSpeedY, params.meshWidth);

    const Real firstTimestep = 1e-5;
    const Real maxTimestep = params.finalTime * 1e-2;
    const Real cflTimestep = std::min(timestepX, timestepY);
    const Real timestep = (time > 0? std::min(maxTimestep, cflTimestep): firstTimestep);


    const CellIndex drillCell = findDrillCell(saturationsWater.rows(), saturationsWater.cols());
    const CellIndex wellCell = findWellCell(saturationsWater.rows(), saturationsWater.cols());

    const Real oldSaturationWellCell = wellCell(saturationsWater);
    const Matrix saturationDivergences = computeSaturationDivergences(fluxFunctionFactors, fluxesX, fluxesY, params.meshWidth);
    saturationsWater -= timestep * saturationDivergences;
    const Real divergenceSaturationWellCell = timestep * wellCell(saturationDivergences);
    const Real outflowWellCell = clamp(wellCell(saturationsWater), 0, 1) * timestep * std::abs(params.inflowPerUnitDepthWater(time)) / (params.porosity * std::pow(params.meshWidth, 2));
    wellCell(saturationsWater) -= outflowWellCell;

    wellCell(saturationsWater) = clamp(wellCell(saturationsWater), 0, 1);
    //wellCell(saturationsWater) += timestep * wellCell(saturationDivergences) * oldSaturationWell;


    drillCell(saturationsWater) = 1;
    time += timestep;
    LOGGER->debug("old sat well = {}", oldSaturationWellCell);
    LOGGER->debug("sat water =\n{}", saturationsWater);

    LOGGER->debug("flux function factors =\n{}", fluxFunctionFactors);
    LOGGER->debug("pressures =\n{}", pressures);

    LOGGER->debug("pressure dx = \n{}", pressureDerivativesX);
    LOGGER->debug("darcy vx =\n{}", darcyVelocitiesX);

    LOGGER->debug("pressure dy = \n{}", pressureDerivativesY);
    LOGGER->debug("darcy vy = \n{}", darcyVelocitiesY);

    LOGGER->debug("inflow = {}", std::abs(params.inflowPerUnitDepthWater(time)));
    LOGGER->debug("meshwidth = {}", params.meshWidth);

    LOGGER->debug("fluxes x =\n{}", fluxesX);
    LOGGER->debug("fluxes y =\n{}", fluxesY);

    //LOGGER->debug("advection velocities x {}", advectionVelocitiesX);

    LOGGER->debug("max advection speed x = {}", maximumAdvectionSpeedX);
    LOGGER->debug("max advection speed y = {}", maximumAdvectionSpeedY);

    LOGGER->debug("CFL timestep = {}", timestep);

    LOGGER->debug("div sat water =\n{}", saturationDivergences);

    LOGGER->debug("timestep = {}", timestep);
    LOGGER->info("Well cell div sat = {}", divergenceSaturationWellCell);
    LOGGER->info("Well cell outflow water = {}", outflowWellCell);
    LOGGER->info("Well cell sat = {}", wellCell(saturationsWater));
    LOGGER->info("South neighbor of well cell sat = {}", wellCell.neighbor(CellIndex::Direction::SOUTH)(saturationsWater));
    LOGGER->info("West neighbor of well cell sat = {}", wellCell.neighbor(CellIndex::Direction::WEST)(saturationsWater));

    if (wellCell(saturationsWater) < 0) {
        throw std::logic_error("well sat < 0!");
    } else if (wellCell(saturationsWater) > 1) {
        throw std::logic_error("well sat > 1!");
    }

}

Matrix clamp(ConstMatrixRef x, const Real minVal, const Real maxVal) {
    return x.unaryExpr([=] (const Real x) -> Real {return clamp(x, minVal, maxVal);});
}