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

CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const int numberOfRows,
      const int numberOfCols) {

    return {fromCell.linearIndex(numberOfRows),
            toCell.linearIndex(numberOfRows)};
}


void adaptPressureGradientsAtWell(const Real wellPressureNow, ConstMatrixRef pressures, MatrixRef pressureDerivativesX, MatrixRef pressureDerivativesY, const Real meshWidth) {
    const CellIndex wellCell = findWellCell(pressures.rows(), pressures.cols());
    const CellIndex eastDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::EAST);
    const CellIndex northDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::NORTH);
    const CellIndex southDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::SOUTH);
    const CellIndex westDerivativeOfWellCell = centerIndexToBorderIndex(wellCell, CellIndex::Direction::WEST);


    eastDerivativeOfWellCell(pressureDerivativesX) = (wellPressureNow - wellCell(pressures)) / meshWidth;
    northDerivativeOfWellCell(pressureDerivativesY) = (wellPressureNow - wellCell(pressures)) / meshWidth;

}



SparseMatrix assemblePressureSystemWithBC(ConstMatrixRef totalMobilities) {
    LOGGER->debug("Starting to assemble system");

    const int numberOfRows = totalMobilities.rows();
    const int numberOfCols = totalMobilities.cols();


    const int numberOfPairs = numberOfCols * numberOfRows;
    SparseMatrix transmissibilities(numberOfPairs, numberOfPairs);


    transmissibilities.reserve(Eigen::VectorXi::Constant(transmissibilities.cols(), 6));

    CellIndex myself = {0, 0};

    for (myself.j = 0; myself.j < numberOfCols; ++myself.j) {
        for (myself.i = 0; myself.i < numberOfRows; ++myself.i) {


            const CellIndex meToMyself = pressureToTransmissibilityIndex(myself, myself, numberOfRows, numberOfCols);

            constexpr static std::array<CellIndex::Direction, 2> directionsToCheck = {
                  CellIndex::Direction::WEST, CellIndex::Direction::NORTH
            };



            for (const CellIndex::Direction direction: directionsToCheck) {


                if (unlikely(!myself.hasNeighbor(direction, numberOfRows, numberOfCols))) {
                    continue;
                }

                const CellIndex neighbor = myself.neighbor(direction);

                const CellIndex neighborToThemselves = pressureToTransmissibilityIndex(neighbor, neighbor, numberOfRows,
                                                                                       numberOfCols);
                const CellIndex meToNeighbor = pressureToTransmissibilityIndex(myself, neighbor, numberOfRows,
                                                                               numberOfCols);
                const CellIndex neighborToMe = meToNeighbor.transpose();

                const Real currentTransmissibility = computeTransmissibility(totalMobilities, myself, neighbor);

                meToMyself(transmissibilities) += currentTransmissibility;

                neighborToThemselves(transmissibilities) += currentTransmissibility;

                meToNeighbor(transmissibilities) -= currentTransmissibility;
                neighborToMe(transmissibilities) -= currentTransmissibility;

            }
        }
    }

    transmissibilities.makeCompressed();

    LOGGER->debug("Assembled system");
    return transmissibilities;
}


void adaptRhsForPressure(const Real sourceAtWellNow, const Real sourceAtDrillNow, VectorRef rhs, const int numberOfRows,
                         const int numberOfCols) {
    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);

    const int drillCellIndex = drillCell.linearIndex(numberOfRows);

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);
    const int wellCellIndex = wellCell.linearIndex(numberOfRows);

    const CellIndex referenceCell = findReferenceCell(numberOfRows, numberOfCols);
    const int referenceCellIndex = referenceCell.linearIndex(numberOfRows);

    rhs(wellCellIndex) = -std::abs(sourceAtWellNow);
    rhs(drillCellIndex) = +std::abs(sourceAtDrillNow);
    rhs(referenceCellIndex) = 0;
}


Vector projectSourcesIntoRange(ConstMatrixRef sources) {

    LOGGER->debug("Starting to project sources into range");
    Vector sourcesProjectedIntoRange(sources.size() + 1);
    sourcesProjectedIntoRange << Eigen::Map<const Vector>(sources.data(), sources.size()), 0;

    sourcesProjectedIntoRange.head(sources.size()).array() -= sourcesProjectedIntoRange.head(sources.size()).array().mean();

    LOGGER->debug("Projected sources into range");

    return sourcesProjectedIntoRange;
}


Vector solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstVectorRef rhs, ConstVectorRef pressureGuess) {
    LOGGER->debug("Solving system");

    Eigen::ConjugateGradient<SparseMatrix, Eigen::Upper | Eigen::Lower> solver;
    solver.setTolerance(1e-6);
    solver.compute(transmissibilities);

    LOGGER->debug("transmissibilities {}", transmissibilities);

    const Vector result = solver.solveWithGuess(rhs, pressureGuess);

    LOGGER->debug("System solved");
    //LOGGER->debug("result = {}", result);

    /*if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solver didn't converge!");
    }*/
    //std::cout << "Estimated error = " << solver.error() << "\n";

    return result;
}


