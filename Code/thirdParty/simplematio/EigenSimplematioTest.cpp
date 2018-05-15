#include <Eigen/Dense>
#include "EigenSimplematio.hpp"
#include <iostream>
#include <string>
#include <complex>
#include <stdexcept>

using namespace SMIO;


template <typename Matrix>
void
test(
    const std::string& testName,
    Matrix& M
) {
    const std::string fileName("Test " + testName + ".mat");


    EigenMatFile matFile(fileName.c_str());
    matFile.writeVariable("M", M);
    std::cout << "Matrix written\n";
    std::cout << "Written Matrix was\n" << M << "\n";


    Matrix readM(matFile.readVariable<Matrix>("M"));
    std::cout << "Matrix read\n";

    std::cout << "Read Matrix was\n" << readM << "\n";


    std::cout << "Test " + testName;
    for (int j(0); j < M.cols(); ++j) {
        for (int i(0); i < M.rows(); ++i) {
            if (M.coeffRef(i, j) != readM.coeffRef(i, j)) {
                std::cout << " failed.\n";
                return;
            }
        }
    }

    std::cout << " passed.\n";

    std::cout.flush();

}

template <typename Scalar>
void
testScalar(
    const std::string testName,
    Scalar& s
) {
    std::cout << "Test " << testName << " ";

    const std::string fileName(testName + ".mat");
    EigenMatFile matFile(fileName.c_str());
    matFile.writeVariable("s", s);
    if (matFile.readVariable<Scalar>("s") == s) {
        std::cout << "passed.\n";
    } else {
        std::cout << "failed.\n";
    }
}


int
main(
) {
    Matrix realM(2, 2);
    realM << 1, 2, 3, 4;
    SparseMatrix sparseM(2, 2);
    sparseM.coeffRef(0, 0) = 1;
    sparseM.coeffRef(0, 1) = 2;
    sparseM.coeffRef(1, 0) = 3;
    sparseM.makeCompressed();

    Vector realV(6);
    realV << 1, 2, 3, 4, 5, 6;

    const ComplexScalar imaginaryUnit(0, 1);
    ComplexMatrix complexM(2, 2);
    complexM << imaginaryUnit, 0, imaginaryUnit, 1;
    ComplexVector complexV(complexM.col(0));

    SparseComplexMatrix sparseComplexM(2, 2);
    sparseComplexM.coeffRef(0, 0) = imaginaryUnit;
    sparseComplexM.coeffRef(0, 1) = -imaginaryUnit;
    sparseComplexM.coeffRef(1, 0) = 2;
    sparseComplexM.makeCompressed();


    test<Matrix>("dense real matrix", realM);
    test<Vector>("dense real vector", realV);
    test<SparseMatrix>("sparse real matrix", sparseM);
    test<ComplexMatrix>("dense complex matrix", complexM);
    test<ComplexVector>("dense complex vector", complexV);
    test<SparseComplexMatrix>("sparse complex matrix", sparseComplexM);
    testScalar<SMIO::Scalar>("real scalar", realV(0));
    testScalar<SMIO::ComplexScalar>("complex scalar", complexV(0));


    EigenMatFile matFile("update test.mat");
    matFile.writeVariable("v", realV);
    matFile.writeVariable("M", realM);
    realV.setZero();
    matFile.writeVariable("v", realV);
    matFile.close();

    std::cout << "Test update matrix ";
    matFile.open("update test.mat");
    if (matFile.readVariable<Vector>("v") == realV && matFile.readVariable<Matrix>("M") == realM) {
        std::cout << "passed.\n";
    } else {
        std::cout << "failed.\n";
    }


    std::cout << "Test dimension checking";
    matFile.open("dimension checking test.mat");
    try {
        Eigen::Matrix<SMIO::Scalar, 3, 3> badDimensions(matFile.readVariable<Eigen::Matrix<SMIO::Scalar, 3, 3>>("v"));
        std::cout << " failed.\n";
    } catch (const std::runtime_error& badDimensionError) {
        std::cout << " passed.\n";
    }
    return 0;
}

