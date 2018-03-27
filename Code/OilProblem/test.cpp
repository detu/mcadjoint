//
// Created by Stefano Weidmann on 24.03.18.
//

#include <acutest.h>
#include "oilproblem.hpp"
#include <iostream>


void testPressurePoisson() {
    const int n = 100;

    const Matrix lambdas = Matrix::Random(n, n).unaryExpr([] (Real x) {return 10*x*x + 1;});
    if (!(lambdas.array() > 0).all()) {
        std::cerr << "Not all lambdas positive!\n";
        TEST_CHECK(false);
    }
    Matrix sources(n, n);
    sources.setZero();
    sources(0, n-1) = 1;

    const SparseMatrix transmissibilities = assembleTransmissibilityMatrix(lambdas);

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
                          computeTransmissibility(lambdas, {i, j}, {i + 1, j}) *
                          (pressures(i, j) - pressures(i + 1, j));
                }

                if (i > 0) {
                    shouldBeSource += computeTransmissibility(lambdas, {i, j}, {i - 1, j}) *
                                      (pressures(i, j) - pressures(i - 1, j));
                }

                if (j < n - 1) {
                    shouldBeSource +=
                          computeTransmissibility(lambdas, {i, j}, {i, j + 1}) *
                          (pressures(i, j) - pressures(i, j + 1));
                }

                if (j > 0) {
                    shouldBeSource +=
                          computeTransmissibility(lambdas, {i, j}, {i, j - 1}) *
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

    std::cout << "error = " << error << "\n";
    TEST_CHECK(error < 1e-5);
}

TEST_LIST = {
      {"solvePressurePoissonTest", testPressurePoisson},
      {nullptr, nullptr}
};

