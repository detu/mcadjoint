//
// Created by Stefano Weidmann on 24.03.18.
//

#include <acutest.h>
#include "oilproblem.hpp"
#include <cmath>
#include <iostream>

void testPressurePoisson() {
    const int n = 10;

    const Matrix lambdas = Matrix::Random(n, n).unaryExpr([] (Real x) {return x+5;});
    Matrix sources(n, n);
    sources.setZero();
    sources(0, n-1) = 1;

    const SparseMatrix transmissibilities = assembleTransmissibilityMatrix(lambdas);



    const Matrix pressures = solvePressurePoissonProblem(transmissibilities, sources);

    Real error = -1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double shouldBeSource = 0;
            if (i < n-1) {
                shouldBeSource +=
                      computeTransmissibility(lambdas, {i, j}, {i + 1, j}) * (pressures(i, j) - pressures(i + 1, j));
            }

            if (i > 0) {
                shouldBeSource += computeTransmissibility(lambdas, {i, j}, {i-1, j}) * (pressures(i, j) - pressures(i-1, j));
            }

            if (j < n-1) {
                shouldBeSource +=
                      computeTransmissibility(lambdas, {i, j}, {i, j + 1}) * (pressures(i, j) - pressures(i, j + 1));
            }

            if (j > 0) {
                shouldBeSource +=
                      computeTransmissibility(lambdas, {i, j}, {i, j - 1}) * (pressures(i, j) - pressures(i, j - 1));
            }

            error = std::max(error, std::abs(sources(i, j) - shouldBeSource));

        }
    }

    std::cout << "error = " << error << "\n";
    TEST_CHECK(error < 1e-5);
}

TEST_LIST = {
      {"solvePressurePoissonTest", testPressurePoisson},
      {nullptr, nullptr}
};

