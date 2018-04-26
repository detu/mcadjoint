//
// Created by Stefano Weidmann on 16.04.18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/dev.hpp>

SparseMatrix computePressureResidualsDerivedByPressure(const SparseMatrix& pressureSystem) {
    return pressureSystem;
}

SparseMatrix computePressureResidualsDerivedBySaturationWater(ConstMatrixRef pressures, ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater) {
    const int numberOfRows = pressures.rows();
    const int numberOfCols = pressures.cols();
    const int numberOfPairs = numberOfCols * numberOfRows;

    SparseMatrix derivatives(numberOfPairs, numberOfPairs);
    derivatives.reserve(Eigen::VectorXi::Constant(derivatives.cols(), 5));

    CellIndex myself = {0, 0};
    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);

    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {

            constexpr std::array<CellIndex::Direction, 4> directionsToCheck = {
                  CellIndex::Direction::SOUTH, CellIndex::Direction::NORTH, CellIndex::Direction::EAST, CellIndex::Direction::WEST
            };

            const bool iAmTheCellInTheWell = myself == wellCell;

            if (iAmTheCellInTheWell) {
                continue;
            }

            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows);
            const Real myDerivativeOfMobility = myself(totalMobilitiesDerivedBySaturationsWater);
            const Real myMobility = myself(totalMobilities);



            for (const auto direction: directionsToCheck) {
                if (!myself.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);
                const Real neighborDerivativeOfMobility = neighbor(totalMobilitiesDerivedBySaturationsWater);
                const Real neighborMobility = neighbor(totalMobilities);


                const Real myPressure = myself(pressures);
                const Real neighborPressure = neighbor(pressures);


                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);


                const Real myMobilityFraction = myMobility / (myMobility + neighborMobility);
                const Real neighborMobilityFraction = neighborMobility / (myMobility + neighborMobility);

                const Real pressureDifference = myPressure - neighborPressure;


                meToMyself(derivatives) += std::pow(neighborMobilityFraction, 2) * pressureDifference;
                meToNeighbor(derivatives) = 2 * neighborDerivativeOfMobility * pressureDifference * std::pow(myMobilityFraction, 2);

            }

            meToMyself(derivatives) *= 2 * myDerivativeOfMobility;


        }
    }

    derivatives.makeCompressed();
    return derivatives;
}


Matrix computeFluxFunctionFactorDerivatives(ConstMatrixRef saturationsWater, const Real porosity, const Real dynamicViscosityWater, const Real dynamicViscosityOil) {
    const auto saturationsWaterArray = saturationsWater.array();
    const auto saturationsOilArray = 1 - saturationsWaterArray;

    return (2 * dynamicViscosityOil * dynamicViscosityWater * saturationsOilArray * saturationsWaterArray).cwiseQuotient(
          porosity * (dynamicViscosityWater * saturationsOilArray.square() + dynamicViscosityOil * saturationsWaterArray.square())).matrix();
}


static inline bool checkWhetherFluxGoesToNeighbor(const CellIndex cell, const CellIndex::Direction direction, ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY) {
    switch (direction) {
        case CellIndex::Direction::EAST: {
            return centerIndexToBorderIndex(cell, direction)(darcyVelocitiesX) > 0;
        }
        case CellIndex::Direction::WEST: {
            return centerIndexToBorderIndex(cell, direction)(darcyVelocitiesX) < 0;
        }
        case CellIndex::Direction::NORTH: {
            return centerIndexToBorderIndex(cell, direction)(darcyVelocitiesY) > 0;
        }
        case CellIndex::Direction::SOUTH: {
            return centerIndexToBorderIndex(cell, direction)(darcyVelocitiesY) < 0;
        }
    }
}


SparseMatrix computeSaturationWaterResidualsDerivedBySaturationWater(ConstMatrixRef fluxFunctionDerivatives,
                                                               ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
                                                               const Real timestep, const Real meshWidth
) {
    const int numberOfRows = fluxFunctionDerivatives.rows();
    const int numberOfCols = fluxFunctionDerivatives.cols();
    const int numberOfPairs = numberOfRows * numberOfCols;

    SparseMatrix derivatives(numberOfPairs, numberOfPairs);

    derivatives.reserve(Eigen::VectorXi::Constant(derivatives.cols(), 5));

    CellIndex myself = {0, 0};

    const Real referenceVelocity = meshWidth / timestep;


    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {

            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows);


            constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
                  CellIndex::Direction::EAST, CellIndex::Direction::WEST, CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH
            };

            for (const auto direction: directionsToCheck) {
                if (!myself.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);
                const bool fluxGoesToNeighbor = checkWhetherFluxGoesToNeighbor(myself, direction, darcyVelocitiesX,
                                                                               darcyVelocitiesY);

                const Real relevantDarcyVelocity = getDerivativeAtCellBorder(myself, darcyVelocitiesX, darcyVelocitiesY, direction);

                if (fluxGoesToNeighbor) {
                    meToMyself(derivatives) += relevantDarcyVelocity;
                } else {
                    meToNeighbor(derivatives) =
                          neighbor(fluxFunctionDerivatives) * relevantDarcyVelocity / referenceVelocity;
                }
            }


            meToMyself(derivatives) *= myself(fluxFunctionDerivatives) / referenceVelocity;
            meToMyself(derivatives) -= 1;

        }
    }

    derivatives.makeCompressed();
    return derivatives;
}

SparseVector computeCostFunctionDerivedByPressure(const Real computedPressureAtDrill, const Real measuredPressureAtDrill, const int numberOfRows, const int numberOfCols) {
    SparseVector derivative(numberOfCols * numberOfRows);

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    derivative.coeffRef(drillCell.linearIndex(numberOfRows)) = 2 * (computedPressureAtDrill - measuredPressureAtDrill);

    return derivative;
}

SparseVector computePressurePartOfDiagonalBlockTimesCostDerivedByState(const SparseMatrix& pressureResidualsByPressure, ConstMatrixRef pressures, const Real measuredPressureAtDrill) {
    const CellIndex drillCell = findDrillCell(pressures.rows(), pressures.cols());
    const Real computedPressureAtDrill = drillCell(pressures);

    const SparseVector costFunctionDerivedByPressure = computeCostFunctionDerivedByPressure(computedPressureAtDrill, measuredPressureAtDrill, pressures.rows(), pressures.cols());
    return pressureResidualsByPressure * costFunctionDerivedByPressure;
}

SparseMatrix computeSaturationWaterResidualsDerivedByPressure(ConstMatrixRef pressureSystem, ConstMatrixRef fluxFunctionFactors,
                                                        ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
                                                        ConstMatrixRef mobilities,
                                                        const Real timestep, const Real meshWidth) {

    const int numberOfRows = fluxFunctionFactors.rows();
    const int numberOfCols = fluxFunctionFactors.cols();
    const int numberOfPairs = numberOfRows * numberOfCols;

    SparseMatrix derivatives(numberOfPairs, numberOfPairs);

    derivatives.reserve(Eigen::VectorXi::Constant(derivatives.cols(), 5));

    CellIndex myself = {0, 0};
    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);

    const Real discretizationFactor = timestep / std::pow(meshWidth, 2);


    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {

            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows);
            const bool iAmTheCellInTheWell = myself == wellCell;


            constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
                  CellIndex::Direction::EAST, CellIndex::Direction::WEST, CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH
            };


            for (const auto direction: directionsToCheck) {
                if (!myself.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);


                const bool fluxGoesToNeighbor = checkWhetherFluxGoesToNeighbor(myself, direction, darcyVelocitiesX,
                                                                               darcyVelocitiesY);
                const Real myFluxFunctionFactor = myself(fluxFunctionFactors);
                const Real neighborFunctionFactor = myself(fluxFunctionFactors);

                const Real upwindFluxFunctionFactor = fluxGoesToNeighbor? myFluxFunctionFactor: neighborFunctionFactor;
                const Real borderTransmissibility = iAmTheCellInTheWell? computeTransmissibility(mobilities, myself, neighbor): -meToNeighbor(pressureSystem);

                const Real borderContribution = borderTransmissibility * upwindFluxFunctionFactor * discretizationFactor;

                meToMyself(derivatives) += borderContribution;
                meToNeighbor(derivatives) = -borderContribution;

            }

        }
    }

    derivatives.makeCompressed();
    return derivatives;

}

inline static Real hmeanDerivedBySecond(const Real a, const Real b) {
    return 2*std::pow(a/(a+b), 2);
}

SparseMatrix computePressureResidualsByLogPermeability(ConstMatrixRef pressures, ConstMatrixRef totalMobilities) {
    const int numberOfRows = pressures.rows();
    const int numberOfCols = pressures.cols();
    const int numberOfPairs = numberOfCols * numberOfRows;

    SparseMatrix derivatives(numberOfPairs, numberOfPairs);

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);

    CellIndex myself = {0, 0};

    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {
            const bool iAmTheCellInTheWell = myself == wellCell;

            if (iAmTheCellInTheWell) {
                continue;
            }

            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows);
            const Real myMobility = myself(totalMobilities);


            constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
                  CellIndex::Direction::EAST, CellIndex::Direction::WEST, CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH
            };

            for (const auto direction: directionsToCheck) {
                if (!myself.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);
                const Real mobilityNeighbor = neighbor(totalMobilities);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);

                const Real pressureDifference = myself(pressures) - neighbor(pressures);

                meToMyself(derivatives) += pressureDifference * hmeanDerivedBySecond(mobilityNeighbor, myMobility);
                meToNeighbor(derivatives) = pressureDifference * hmeanDerivedBySecond(myMobility, mobilityNeighbor) * mobilityNeighbor;

            }

            meToMyself(derivatives) *= myMobility;
        }
    }

    return derivatives;

}

SparseMatrix computeSaturationsWaterResidualsByLogPermeability(ConstMatrixRef fluxesX, ConstMatrixRef fluxesY, ConstMatrixRef mobilities, const Real timestep, const Real meshWidth) {
    const int numberOfRows = mobilities.rows();
    const int numberOfCols = mobilities.cols();
    const int numberOfPairs = numberOfCols * numberOfRows;

    SparseMatrix derivatives(numberOfPairs, numberOfPairs);

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);

    CellIndex myself = {0, 0};

    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {
            const bool iAmTheCellInTheWell = myself == wellCell;

            if (iAmTheCellInTheWell) {
                continue;
            }

            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows);
            const Real myMobility = myself(mobilities);


            constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
                  CellIndex::Direction::EAST, CellIndex::Direction::WEST, CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH
            };

            constexpr static std::array<int, 4> signsForDirections = {
                  +1, -1, +1, -1
            };

            for (int directionIndex = 0; directionIndex < 4; ++directionIndex) {
                const auto direction = directionsToCheck[directionIndex];
                const int directionSign = signsForDirections[directionIndex];
                if (!myself.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const Real signedFlux = directionSign * getDerivativeAtCellBorder(myself, fluxesX, fluxesY, direction);

                const CellIndex neighbor = myself.neighbor(direction);
                const Real mobilityNeighbor = neighbor(mobilities);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);


                meToMyself(derivatives) += signedFlux * hmeanDerivedBySecond(mobilityNeighbor, myMobility);
                meToNeighbor(derivatives) = signedFlux * hmeanDerivedBySecond(myMobility, mobilityNeighbor) * timestep / meshWidth;

            }

            meToMyself(derivatives) *= timestep / meshWidth;
        }
    }

    return derivatives;

}