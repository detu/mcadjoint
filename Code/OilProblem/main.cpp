//
// Created by Stefano Weidmann on 24.03.18.
//

#include "oilproblem.hpp"
#include <iostream>

int main() {
    std::cout << "Test transmission matrix\n------------------------\n";
    const int n = 5;
    std::cout << "n = " << n << "\n";

    Matrix lambdas(n, n);

    lambdas.setRandom();
    lambdas.array() += 5;

    std::cout << "lambdas = \n" << lambdas << "\n";

    const Matrix transmissibilities = assembleTransmissibilityMatrix(lambdas);

    Eigen::JacobiSVD<Matrix> svd(transmissibilities);
    double cond = svd.singularValues()(0)
                  / svd.singularValues()(svd.singularValues().size()-1);

    std::cout << "cond = " << cond << "\n";
    std::cout << "transmissibilities = \n" << transmissibilities;

    return 0;
}
