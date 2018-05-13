//
// Created by Stefano Weidmann on 27.03.18.
//

#include "logging.hpp"
#include <stefCommonHeaders/assert.h>
#include "specialCells.hpp"
#include "pressure.hpp"

static inline Real harmonicMean(const Real a, const Real b) {
    return 2 * a * b / (a + b);
}

Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell) {
    // transmissibility is harmonic mean of total mobility coefficients

    const Real lambdaFrom = fromCell(totalMobilities);
    const Real lambdaTo = toCell(totalMobilities);

    return harmonicMean(lambdaFrom, lambdaTo);
}


Real computeTransmissibility(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const CellIndex& fromCell, const CellIndex& toCell) {
    const Real lambdaFrom = computeTotalMobility(dynamicViscosityOil, dynamicViscosityWater, permeabilities, saturationsWater, fromCell);
    const Real lambdaTo = computeTotalMobility(dynamicViscosityOil, dynamicViscosityWater, permeabilities, saturationsWater, toCell);

    return harmonicMean(fromCell, toCell);

}

SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities) {
    //log()->debug("Starting to assemble system");

    const int numberOfRows = totalMobilities.rows();
    const int numberOfCols = totalMobilities.cols();


    const int numberOfPairs = numberOfCols * numberOfRows;
    SparseMatrix transmissibilities(numberOfPairs, numberOfPairs);

    //const CellIndex referenceCell = findReferenceCell(numberOfRows, numberOfCols);

    transmissibilities.reserve(Eigen::VectorXi::Constant(transmissibilities.cols(), 5));

    CellIndex myself = {0, 0};
    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const CellIndex cellWithTheKnownPressure = wellCell;

    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {

            const bool iAmTheCellWithTheKnownPressure = myself == cellWithTheKnownPressure;




            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows);

            if (iAmTheCellWithTheKnownPressure) {
                meToMyself(transmissibilities) = 1;
            }

            constexpr static std::array<CellIndex::Direction, 2> directionsToCheck = {
                  CellIndex::Direction::EAST, CellIndex::Direction::SOUTH
            };



            for (const CellIndex::Direction direction: directionsToCheck) {

                if (!myself.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);

                const CellIndex neighborToThemselves = pressureToTransmissibilityIndex(neighbor, neighbor, numberOfRows);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows);
                const CellIndex neighborToMe = meToNeighbor.transpose();

                const Real currentTransmissibility = computeTransmissibility(totalMobilities, myself, neighbor);

                if (!iAmTheCellWithTheKnownPressure) {
                    meToMyself(transmissibilities) += currentTransmissibility;
                    meToNeighbor(transmissibilities) -= currentTransmissibility;
                }

                if (neighbor != cellWithTheKnownPressure) {
                    neighborToThemselves(transmissibilities) += currentTransmissibility;
                    neighborToMe(transmissibilities) -= currentTransmissibility;
                }

            }
        }
    }

    transmissibilities.makeCompressed();

    //log()->debug("Assembled system");
    return transmissibilities;
}


Vector computeRhsForPressureSystem(const Real sourceAtDrillNow, const int numberOfRows, const int numberOfCols) {
    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);

    const int drillCellIndex = drillCell.linearIndex(numberOfRows);


    Vector rhs(Vector::Zero(numberOfRows*numberOfCols));
    rhs(drillCellIndex) = +std::abs(sourceAtDrillNow);

    return rhs;
}





Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities,
                                   ConstVectorRef rhs) {
    log()->debug("Solving system");

    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(transmissibilities);

    log()->debug("transmissibilities =\n{}", transmissibilities);
    log()->debug("rhs =\n{}", Vector(rhs));

    const Vector result = solver.solve(rhs);


    log()->debug("System solved");
    log()->debug("result = {}", result);

    /*if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solver didn't converge!");
    }*/
    //std::cout << "Estimated error = " << solver.error() << "\n";

    return result;
}

Matrix computeTotalMobilities(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater) {
    return permeabilities.array() * (saturationsWater.array().square() / dynamicViscosityWater + (1.0 - saturationsWater.array()).square() / dynamicViscosityOil);
}

Real computeTotalMobility(const Real dynamicViscosityOil, const Real dynamicViscosityWater, ConstMatrixRef permeabilities, ConstMatrixRef saturationsWater, const CellIndex entry) {
    const Real saturation = entry(saturationsWater);
    return entry(permeabilities) * (std::pow(saturation, 2) / dynamicViscosityWater + std::pow(1 - saturation, 2) / dynamicViscosityOil);
}

