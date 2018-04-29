//
// Created by Stefano Weidmann on 27.03.18.
//

#include "oilProblem.hpp"
#include "logging.hpp"
#include <stefCommonHeaders/assert.h>
#include "specialCells.hpp"

#ifdef __GNUC__
    #define likely(x)       __builtin_expect((x),1)
    #define unlikely(x)     __builtin_expect((x),0)
#else
    #define likely(x) x
    #define unlikely(x) x
#endif



Real computeTransmissibility(ConstMatrixRef totalMobilities, const CellIndex& fromCell, const CellIndex& toCell) {
    // transmissibility is harmonic mean of total mobility coefficients

    const double lambdaFrom = fromCell(totalMobilities);
    const double lambdaTo = toCell(totalMobilities);

    return 2.0 * lambdaFrom * lambdaTo / (lambdaFrom + lambdaTo);
}

CellIndex pressureToTransmissibilityIndex(const CellIndex& fromCell, const CellIndex& toCell, const int numberOfRows) {

    return {fromCell.linearIndex(numberOfRows),
            toCell.linearIndex(numberOfRows)};
}


void adaptPressureGradientsAtWell(const Real inflowNow, ConstMatrixRef mobilities, ConstMatrixRef pressures, MatrixRef pressureDerivativesX, MatrixRef pressureDerivativesY, const Real meshWidth) {
    const CellIndex wellCell = findWellCell(pressures.rows(), pressures.cols());
    const CellIndex eastDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::EAST);
    const CellIndex northDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::NORTH);
    const CellIndex southDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::SOUTH);
    const CellIndex westDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::WEST);


    eastDerivativeOfWellCell(pressureDerivativesX) = -inflowNow / (4 * wellCell(mobilities) * meshWidth);
    northDerivativeOfWellCell(pressureDerivativesY) = eastDerivativeOfWellCell(pressureDerivativesX);

}



SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities) {
    //LOGGER->debug("Starting to assemble system");

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

                if (unlikely(!myself.hasNeighbor(direction, numberOfRows, numberOfCols))) {
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

    //LOGGER->debug("Assembled system");
    return transmissibilities;
}


SparseVector computeRhsForPressureSystem(const Real sourceAtDrillNow, const int numberOfRows, const int numberOfCols) {
    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);

    const int drillCellIndex = drillCell.linearIndex(numberOfRows);


    SparseVector rhs(numberOfRows*numberOfCols);
    rhs.coeffRef(drillCellIndex) = +std::abs(sourceAtDrillNow);

    return rhs;
}


Vector projectSourcesIntoRange(ConstMatrixRef sources) {

    LOGGER->debug("Starting to project sources into range");
    Vector sourcesProjectedIntoRange(sources.size() + 1);
    sourcesProjectedIntoRange << Eigen::Map<const Vector>(sources.data(), sources.size()), 0;

    sourcesProjectedIntoRange.head(sources.size()).array() -= sourcesProjectedIntoRange.head(sources.size()).array().mean();

    LOGGER->debug("Projected sources into range");

    return sourcesProjectedIntoRange;
}


Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities, const SparseVector& rhs) {
    LOGGER->debug("Solving system");

    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(transmissibilities);

    LOGGER->debug("transmissibilities =\n{}", transmissibilities);
    LOGGER->debug("rhs =\n{}", rhs);

    const Vector result = solver.solve(rhs);

    LOGGER->debug("System solved");
    //LOGGER->debug("result = {}", result);

    /*if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solver didn't converge!");
    }*/
    //std::cout << "Estimated error = " << solver.error() << "\n";

    return result;
}



