//
// Created by Stefano Weidmann on 16.04.18.
//

#include "derivativesForAdjoint.hpp"
#include "darcyVelocity.hpp"
#include "specialCells.hpp"
#include "pressure.hpp"
#include <stefCommonHeaders/dev.hpp>

inline static Real hmeanDerivedBySecond(const Real a, const Real b) {
    return 2*std::pow(a/(a+b), 2);
}


SparseMatrix computePressureResidualsDerivedByPressure(const SparseMatrix& pressureSystem) {
    return pressureSystem;
}

Matrix computeTotalMobilitiesDerivedBySaturationsWater(ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const Real dynamicViscosityOil, const Real dynamicViscosityWater) {
    const int numberOfRows = permeabilities.rows();
    const int numberOfCols = permeabilities.cols();


    const auto saturationsArray = saturationsWater.array();
    const auto permeabilitiesArray = permeabilities.array();

    return (2*permeabilitiesArray * ((saturationsArray - 1) / dynamicViscosityOil + saturationsArray / dynamicViscosityWater)).matrix();
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
          porosity * (dynamicViscosityWater * saturationsOilArray.square() + dynamicViscosityOil * saturationsWaterArray.square()).square()).matrix();
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


SparseMatrix computeSaturationWaterResidualsDerivedBySaturationWater(
      ConstMatrixRef fluxFunctionFactors, ConstMatrixRef fluxFunctionDerivatives,
      ConstMatrixRef darcyVelocitiesX, ConstMatrixRef darcyVelocitiesY,
      ConstMatrixRef pressureDerivativesX, ConstMatrixRef pressureDerivativesY,
      ConstMatrixRef totalMobilities, ConstMatrixRef totalMobilitiesDerivedBySaturationsWater,
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

            constexpr static std::array<int, 4> signs = {+1, -1, +1, -1};

            int signIndex = -1;
            for (const auto direction: directionsToCheck) {
                ++signIndex;
                if (!myself.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const int sign = signs[signIndex];

                const CellIndex neighbor = myself.neighbor(direction);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);
                const bool fluxGoesToNeighbor = checkWhetherFluxGoesToNeighbor(myself, direction, darcyVelocitiesX,
                                                                               darcyVelocitiesY);

                const Real relevantDarcyVelocity = getDerivativeAtCellBorder(myself, darcyVelocitiesX, darcyVelocitiesY, direction);
                const Real upstreamFluxFunction = fluxGoesToNeighbor? myself(fluxFunctionFactors): neighbor(fluxFunctionFactors);
                const Real pressureGradientAtBoundary = getDerivativeAtCellBorder(myself, pressureDerivativesX, pressureDerivativesY, direction);

                const Real myMobility = myself(totalMobilities);
                const Real neighborMobility = neighbor(totalMobilities);

                const Real myMobilityDerivative = myself(totalMobilitiesDerivedBySaturationsWater);
                const Real neighborMobilityDerivative = neighbor(totalMobilitiesDerivedBySaturationsWater);

                const Real mobilityTermMyself = hmeanDerivedBySecond(neighborMobility, myMobility) * myMobilityDerivative;
                const Real mobilityTermNeighbor = hmeanDerivedBySecond(myMobility, neighborMobility) * neighborMobilityDerivative;

                if (fluxGoesToNeighbor) {
                    meToMyself(derivatives) += sign * relevantDarcyVelocity;
                } else {
                    meToNeighbor(derivatives) = neighbor(fluxFunctionDerivatives) * relevantDarcyVelocity;
                }

                meToNeighbor(derivatives) -= mobilityTermNeighbor * pressureGradientAtBoundary;
                meToNeighbor(derivatives) *= sign / referenceVelocity;

                meToMyself(derivatives) -= sign * pressureGradientAtBoundary * mobilityTermMyself;
            }


            meToMyself(derivatives) *= myself(fluxFunctionDerivatives) / referenceVelocity;
            meToMyself(derivatives) -= 1;

        }
    }

    derivatives.makeCompressed();
    return derivatives;
}

Real computeCostFunctionDerivedByPressureEntry(const Real computedPressureAtDrill, const Real measuredPressureAtDrill,
                                               const int numberOfRows, const int numberOfCols,
                                               const CellIndex& derivedByCell) {

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    if (drillCell == derivedByCell) {
        return 2 * (computedPressureAtDrill - measuredPressureAtDrill);
    } else {
        return 0;
    }

}

Real computeCostFunctionDerivedBySaturationsWaterEntry(const int numberOfRows, const int numberOfCols,
                                                       const CellIndex& derivedByCell) {
    return 0;
}



SparseMatrix computeSaturationWaterResidualsDerivedByPressure(const SparseMatrix& pressureSystem, ConstMatrixRef fluxFunctionFactors,
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
                const Real neighborFunctionFactor = neighbor(fluxFunctionFactors);

                const Real upwindFluxFunctionFactor = fluxGoesToNeighbor? myFluxFunctionFactor: neighborFunctionFactor;
                const Real borderTransmissibility = iAmTheCellInTheWell? computeTransmissibility(mobilities, myself, neighbor): -meToNeighbor(pressureSystem);

                const Real borderContribution = borderTransmissibility * upwindFluxFunctionFactor * discretizationFactor;

                meToMyself(derivatives) -= borderContribution;
                meToNeighbor(derivatives) = borderContribution;

            }

        }
    }

    derivatives.makeCompressed();
    return derivatives;

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

SparseMatrix computeSaturationsWaterResidualsByLogPermeability(ConstMatrixRef pressureGradientsX, ConstMatrixRef pressureGradientsY, ConstMatrixRef mobilities, const Real timestep, const Real meshWidth) {
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

                const Real signedPressureGradient = -directionSign * getDerivativeAtCellBorder(myself, pressureGradientsX, pressureGradientsY, direction);

                const CellIndex neighbor = myself.neighbor(direction);
                const Real mobilityNeighbor = neighbor(mobilities);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);


                meToMyself(derivatives) += signedPressureGradient * hmeanDerivedBySecond(mobilityNeighbor, myMobility);
                meToNeighbor(derivatives) = signedPressureGradient * hmeanDerivedBySecond(myMobility, mobilityNeighbor) * timestep / meshWidth * mobilityNeighbor;

            }

            meToMyself(derivatives) *= timestep / meshWidth * myMobility;
        }
    }

    return derivatives;

}

