//
// Created by Stefano Weidmann on 24.03.18.
//

#include <acutest.h>
#include <iostream>
#include "typedefs.hpp"
#include "cellindex.hpp"
#include "specialCells.hpp"
#include "pressure.hpp"
#include "darcyVelocity.hpp"
#include "dumpToMatFile.hpp"

const char* MAT_FILE_NAME = "test.mat";

void testPressurePoisson() {
    const int n = 100;

    const Matrix totalMobilities = Matrix::Random(n, n).unaryExpr([] (Real x) {return 10 * std::abs(x) + 1;});
    if (!(totalMobilities.array() > 0).all()) {
        std::cerr << "Not all totalMobilities positive!\n";
        TEST_CHECK(false);
    }

    const CellIndex wellCell = findWellCell(n, n);
    const CellIndex drillCell = findDrillCell(n, n);

    Matrix sources(Matrix::Zero(n, n));
    wellCell(sources) = -1;

    drillCell(sources) = -wellCell(sources);

    Matrix pressures(n, n);
    pressures.setZero();


    const SparseMatrix transmissibilities = assemblePressureSystemWithBC(totalMobilities);

    Eigen::Map<Vector> pressuresAsVector(pressures.data(), pressures.size());



    const SparseVector rhs = computeRhsForPressureSystem(-1, n, n);

    pressuresAsVector = std::move(solvePressurePoissonProblem(transmissibilities, rhs));
    #ifdef VERBOSE_TESTS
    std::cout << "pressures = " << pressures << "\n";
    std::cout << "pressures as vector = " << pressuresAsVector << "\n";

    std::cout << "Computed sources = " << (transmissibilities * pressuresAsVector).transpose() << "\n";
    std::cout << "Real Sources = " << sources << "\n";
    #endif

    Real error = -1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double shouldBeNegativeSource = 0;

            if (i < n - 1) {
                shouldBeNegativeSource +=
                      computeTransmissibility(totalMobilities, {i, j}, {i + 1, j}) *
                      (pressures(i, j) - pressures(i + 1, j));
            }

            if (i > 0) {
                shouldBeNegativeSource += computeTransmissibility(totalMobilities, {i, j}, {i - 1, j}) *
                                  (pressures(i, j) - pressures(i - 1, j));
            }

            if (j < n - 1) {
                shouldBeNegativeSource +=
                      computeTransmissibility(totalMobilities, {i, j}, {i, j + 1}) *
                      (pressures(i, j) - pressures(i, j + 1));
            }

            if (j > 0) {
                shouldBeNegativeSource +=
                      computeTransmissibility(totalMobilities, {i, j}, {i, j - 1}) *
                      (pressures(i, j) - pressures(i, j - 1));
            }

            const double newError = std::abs(sources(i, j) + shouldBeNegativeSource);
            if (newError > error) {
                #ifdef VERBOSE_TESTS
                std::cout << "New max error at (" << i << ", " << j << ") of " << newError << "\n"
                          << "is " << -shouldBeNegativeSource << " but should be " << sources(i, j) << "\n";
                #endif
                error = newError;
            }

        }
    }

    //std::cout << "error = " << error << "\n";
    TEST_CHECK(error < 1e-3);
}


void testDerivatives() {
    const Real meshWidth = 1e-2;
    const int n = 5;

    const auto f = [] (const Point& p) {
        return std::sin(p.x) * std::cos(p.y);
    };

    const auto dfdx = [] (const Point& p) {
        return std::cos(p.x) * std::cos(p.y);
    };

    const auto dfdy = [] (const Point& p) {
        return -std::sin(p.x) * std::sin(p.y);
    };



    Matrix field(n, n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            field(i, j) = f(CellIndex({i, j}).toCenterPoint(meshWidth));
        }
    }

    const Matrix computedXDerivatives = computeXDerivative(field, meshWidth);
    const Matrix computedYDerivatives = computeYDerivative(field, meshWidth);


    #ifdef VERBOSE_TESTS
    std::cout << "f =\n" << field << "\n\ndx =\n" << computedXDerivatives << "\n\ndy =\n" << computedYDerivatives << "\n";
    #endif



    constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
          CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH, CellIndex::Direction::EAST, CellIndex::Direction::WEST
    };

    CellIndex cell = {0, 0};
    CellIndex greatestErrorLocation;
    CellIndex::Direction greatestErrorDirection;
    Real greatestError = -1;
    Real greatestErrorExactDerivative = NAN;
    Real greatestErrorComputedDerivative = NAN;
    for (cell.j = 0; cell.j < n; ++cell.j) {
        for (cell.i = 0; cell.i < n; ++cell.i) {
            for (const CellIndex::Direction direction: directionsToCheck) {

                const Real computedDerivative = getDerivativeAtCellBorder(cell, computedXDerivatives, computedYDerivatives, direction);
                if (cell.j == 0 && direction == CellIndex::Direction::WEST ||
                    cell.j == n-1 && direction == CellIndex::Direction::EAST ||
                    cell.i == 0 && direction == CellIndex::Direction::NORTH  ||
                    cell.i == n-1 && direction == CellIndex::Direction::SOUTH) {
                    TEST_CHECK_(computedDerivative == 0, "derivative at border not zero but %lf", computedDerivative);
                    continue;
                }


                Real exactDerivative = NAN;
                const Point locationOfBorder = cell.toBorderPoint(meshWidth, direction);
                switch (direction) {
                    case CellIndex::Direction::NORTH:
                    case CellIndex::Direction::SOUTH: {
                        exactDerivative = dfdy(locationOfBorder);
                        break;
                    }

                    case CellIndex::Direction::EAST:
                    case CellIndex::Direction::WEST: {
                        exactDerivative = dfdx(locationOfBorder);
                        break;
                    }
                }
                const Real error = std::abs(exactDerivative - computedDerivative);

                if (error > greatestError) {
                    greatestError = error;
                    greatestErrorLocation = cell;
                    greatestErrorDirection = direction;
                    greatestErrorExactDerivative = exactDerivative;
                    greatestErrorComputedDerivative = computedDerivative;
                }
            }
        }
    }

    TEST_CHECK_(greatestError < 1e-5, "Greatest error %lf is at (%d, %d) in direction %s\nThere the exact derivative is %lf, but computed was %lf",
                greatestError, greatestErrorLocation.i, greatestErrorLocation.j, CellIndex::directionToString(greatestErrorDirection),
                greatestErrorExactDerivative, greatestErrorComputedDerivative);
}

TEST_LIST = {
      {"solvePressurePoissonTest", testPressurePoisson},
      {"derivativesTest", testDerivatives},
      {nullptr, nullptr}
};

