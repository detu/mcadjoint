//
// Created by Stefano Weidmann on 24.03.18.
//

#include "oilProblem.hpp"
#include <iostream>

int main() {
    std::cout << "Test transmission matrix\n------------------------\n";
    const int n = 2;
    std::cout << "n = " << n << "\n";

    Matrix totalMobilities(n, n);

    totalMobilities.setRandom();
    totalMobilities.array() += 5;

    if (n <= 10) {
        std::cout << "totalMobilities = \n" << totalMobilities << "\n";
    }

    Matrix sources;
    sources.resizeLike(totalMobilities);
    sources.setZero();

    const Matrix transmissibilities = assemblePressureSystemWithBC(totalMobilities);

    Eigen::JacobiSVD<Matrix> svd(transmissibilities);
    double cond = svd.singularValues()(0)
                  / svd.singularValues()(svd.singularValues().size()-1);

    std::cout << "cond = " << cond << "\n";
    std::cout << "eigenvalues = " << (transmissibilities.transpose() * transmissibilities).eigenvalues() << "\n";

    if (n <= 10) {
        std::cout << "transmissibilities = \n" << transmissibilities;
    } else {
        std::cout << "too big, don't display transmissibilities\n";
    }
    return 0;
}

