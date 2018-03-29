//
// Created by Stefano Weidmann on 24.03.18.
//

#include <acutest.h>
#include "oilProblem.hpp"
#include <iostream>


void testPressurePoisson() {
    const int n = 100;

    const Matrix totalMobilities = Matrix::Random(n, n).unaryExpr([] (Real x) {return 10 * std::abs(x) + 1;});
    if (!(totalMobilities.array() > 0).all()) {
        std::cerr << "Not all totalMobilities positive!\n";
        TEST_CHECK(false);
    }
    Matrix sources(Matrix::Random(n, n));

    const SparseMatrix transmissibilities = assembleTransmissibilityMatrix(totalMobilities);

    VectorToBeMappedAsMatrix vecmap = solvePressurePoissonProblem(transmissibilities, sources);
    const Matrix& pressures = vecmap.map;

    const CellIndex wellCell = {0, n-1};
    Real error = -1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double shouldBeSource = 0;


            if (CellIndex({i, j}) == wellCell) {
                shouldBeSource = pressures(i, j);
            } else {
                if (i < n - 1) {
                    shouldBeSource +=
                          computeTransmissibility(totalMobilities, {i, j}, {i + 1, j}) *
                          (pressures(i, j) - pressures(i + 1, j));
                }

                if (i > 0) {
                    shouldBeSource += computeTransmissibility(totalMobilities, {i, j}, {i - 1, j}) *
                                      (pressures(i, j) - pressures(i - 1, j));
                }

                if (j < n - 1) {
                    shouldBeSource +=
                          computeTransmissibility(totalMobilities, {i, j}, {i, j + 1}) *
                          (pressures(i, j) - pressures(i, j + 1));
                }

                if (j > 0) {
                    shouldBeSource +=
                          computeTransmissibility(totalMobilities, {i, j}, {i, j - 1}) *
                          (pressures(i, j) - pressures(i, j - 1));
                }
            }

            const double newError = std::abs(sources(i, j) - shouldBeSource);
            if (newError > error) {
                #ifdef VERBOSE_TESTS
                std::cout << "New max error at (" << i << ", " << j << ") of " << newError << "\n"
                          << "is " << shouldBeSource << " but should be " << sources(i, j) << "\n";
                #endif
                error = newError;
            }

        }
    }

    //std::cout << "error = " << error << "\n";
    TEST_CHECK(error < 1e-5);
}


void testDerivatives() {
    const Real meshWidth = 1e-6;
    const int n = 100;

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
            field(i, j) = f(CellIndex({i, j}).toCenterPoint(n, n, meshWidth));
        }
    }

    const Matrix computedXDerivatives = computeXDerivative(field, meshWidth);
    const Matrix computedYDerivatives = computeYDerivative(field, meshWidth);

    constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
          CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH, CellIndex::Direction::EAST, CellIndex::Direction::WEST
    };

    CellIndex cell = {0, 0};
    CellIndex greatestErrorLocation;
    Real greatestError = -1;
    for (cell.j = 0; cell.j < n; ++cell.j) {
        for (cell.i = 0; cell.i < n; ++cell.i) {
            for (const CellIndex::Direction direction: directionsToCheck) {

                const Real computedDerivative = getDerivativeAtCellBorder(cell, computedXDerivatives, computedYDerivatives, direction);
                if (cell.j == 0 && direction == CellIndex::Direction::EAST ||
                    cell.j == n-1 && direction == CellIndex::Direction::WEST ||
                    cell.i == 0 && direction == CellIndex::Direction::NORTH  ||
                    cell.i == n-1 && direction == CellIndex::Direction::SOUTH) {
                    TEST_CHECK_(computedDerivative == 0, "derivative at border not zero but %lf", computedDerivative);
                    continue;
                }


                Real exactDerivative = NAN;
                const Point locationOfBorder = cell.toBorderPoint(n, n, meshWidth, direction);
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
                }
            }
        }
    }

    TEST_CHECK_(greatestError < 1e-5, "Greatest error %lf is at (%d, %d)", greatestError, greatestErrorLocation.i, greatestErrorLocation.j);
}

TEST_LIST = {
      {"solvePressurePoissonTest", testPressurePoisson},
      {"derivativesTest", testDerivatives},
      {nullptr, nullptr}
};

