//
// Created by Stefano Weidmann on 16.05.18.
//

#include "preconditioning.hpp"
#include "pressure.hpp"
#include "logging.hpp"
#include "utils.hpp"
#include "preconditioningOptions.hpp"
#include "specialCells.hpp"

static SparseMatrix computeTargetPressureByPressureMatrix(const SparseMatrix& pressuresByPressures,
                                                          const SparseMatrix& saturationsByPressures,
                                                          const int numberOfRows, const int numberOfCols) {



    SparseMatrixRowMajor target(pressuresByPressures.rows(), pressuresByPressures.cols());
    target.setIdentity();


    return target;
}


void preconditionMatrices(const int numberOfRows, const int numberOfCols,
      SparseMatrix& pressuresByPressures, SparseMatrix& saturationsByPressures,
      SparseMatrix& pressuresBySaturations, SparseMatrix& saturationsBySaturations,
      VectorRef b, const PressureSolver& pressureSolver, const WhichPreconditioner whichPreconditioner) {

    ASSERT(allFinite(b));

    const int systemSize = b.size()/2;
    VectorRef pressurePartOfB = b.head(systemSize);

    switch (whichPreconditioner) {
        default: {
            std::abort();
        }

        case WhichPreconditioner::NONE: {
            log()->info("Preconditioner: None");
            break;
        }

        case WhichPreconditioner::DIAGONAL: {
            log()->info("Preconditioner: inverse diagonal from pressure by pressure");
            const SparseMatrix inversePressureByPressureDiagonal = SparseMatrix(extractInverseDiagonalMatrix(pressuresByPressures));
            pressuresByPressures = pressuresByPressures * inversePressureByPressureDiagonal;
            saturationsByPressures = saturationsByPressures * inversePressureByPressureDiagonal;
            pressurePartOfB = (inversePressureByPressureDiagonal * pressurePartOfB).eval();

            break;

        }


        case WhichPreconditioner::PRESSURE_BY_PRESSURE: {

            const SparseMatrix targetPressureByPressureMatrix = computeTargetPressureByPressureMatrix(pressuresByPressures, saturationsByPressures, numberOfRows, numberOfCols);

            log()->info("Preconditioner: full inverse of  pressure by pressure");
            SparseMatrix saturationsWaterResidualsByPressuresTransposed = saturationsByPressures.transpose();
            log()->debug("nonzeros in saturationsByPressures = {}", saturationsByPressures.nonZeros());
            pressuresByPressures.setZero();
            pressuresByPressures = targetPressureByPressureMatrix;


            const SparseMatrix  saturationsByPressuresPreconditioned = targetPressureByPressureMatrix * pressureSolver.solve(saturationsWaterResidualsByPressuresTransposed).transpose();
            ASSERT(pressureSolver.info() == Eigen::Success);

            saturationsByPressures = saturationsByPressuresPreconditioned;


            saturationsByPressures.makeCompressed();
            ASSERT(saturationsByPressures.valuePtr() != nullptr);
            pressurePartOfB = targetPressureByPressureMatrix * pressureSolver.solve(pressurePartOfB);

            break;
        }

        case WhichPreconditioner::Q_TRANSPOSE_FROM_QR_AND_DIAGONAL: {
            const int numberOfCells = pressuresByPressures.rows();
            log()->info("Preconditioner: Q^T from the QR decomposition of pressure by pressure + diagonal");

            const Eigen::SparseQR<SparseMatrix, Eigen::NaturalOrdering<int>> qrSolver(pressuresByPressures);
            const Eigen::SparseMatrix<Real, Eigen::RowMajor> sortedRMatrix = qrSolver.matrixR();
            pressuresByPressures = SparseMatrix(sortedRMatrix.transpose());

            SparseMatrix inverseDiagonal = SparseMatrix(extractInverseDiagonalMatrix(pressuresByPressures));
            ASSERT(allFinite(inverseDiagonal));
            if (false) {
                for (int colIndex = 0; colIndex < inverseDiagonal.cols(); ++colIndex) {
                    inverseDiagonal.coeffRef(colIndex, colIndex) /= std::max(1.0,
                                                                             pressuresByPressures.col(
                                                                                   colIndex).cwiseAbs().sum());
                }
            }
            pressuresByPressures = (pressuresByPressures * inverseDiagonal).eval();



            SparseMatrix saturationsWaterResidualsByPressuresTransposed = saturationsByPressures.transpose();
            SparseMatrix result(saturationsWaterResidualsByPressuresTransposed.rows(), saturationsWaterResidualsByPressuresTransposed.cols());
            saturationsWaterResidualsByPressuresTransposed.makeCompressed();

            Vector denseCol(numberOfCells);
            for (int colIndex = 0; colIndex < saturationsWaterResidualsByPressuresTransposed.cols(); ++colIndex) {
                denseCol = saturationsWaterResidualsByPressuresTransposed.col(colIndex);
                denseCol = (inverseDiagonal * (qrSolver.matrixQ().transpose() * denseCol)).eval();
                Eigen::SparseVector<Real> sparseCol(denseCol.size());
                for (int elementIndex = 0; elementIndex < denseCol.size(); ++elementIndex) {
                    if (std::abs(denseCol(elementIndex)) > 1e-14) {
                        sparseCol.coeffRef(elementIndex) = denseCol(elementIndex);
                    }
                }
                result.col(colIndex) = sparseCol;
            }
            pressurePartOfB = (inverseDiagonal * (qrSolver.matrixQ().transpose() * pressurePartOfB)).eval();


            saturationsByPressures = result.transpose();
            saturationsByPressures.makeCompressed();
            break;
        }

        case WhichPreconditioner::CHOLESKY_L_INV: {
            log()->info("Preconditioner: Cholesky L^-1");
            const Eigen::SimplicialLLT<SparseMatrix> choleskySolver(pressuresByPressures);

            pressuresByPressures = choleskySolver.matrixL();
            const SparseMatrix inverseDiagonal = SparseMatrix(extractInverseDiagonalMatrix(pressuresByPressures));
            pressuresByPressures = (pressuresByPressures * inverseDiagonal).eval();



            SparseMatrix saturationsWaterResidualsByPressuresTransposed = saturationsByPressures.transpose();
            saturationsWaterResidualsByPressuresTransposed.makeCompressed();

            choleskySolver.matrixL().solveInPlace(saturationsWaterResidualsByPressuresTransposed);
            saturationsByPressures = (inverseDiagonal * saturationsByPressures).eval();

            saturationsByPressures = saturationsWaterResidualsByPressuresTransposed.transpose();


            choleskySolver.matrixL().solveInPlace(pressurePartOfB);
            pressurePartOfB = (inverseDiagonal * pressurePartOfB).eval();
            break;
        }
    }

    ASSERT(allFinite(b));


    if (whichPreconditioner != WhichPreconditioner::NONE) {
        // estimate the variance growth factor gamma_i |b_i| / B + gamma_i^2 on page 6196
        const std::array<const std::array<SparseMatrix*, 2>, 2> jordanBlocksInConsideration = {{
                                                                                                     {&pressuresByPressures, &saturationsByPressures},
                                                                                                     {&pressuresBySaturations, &saturationsBySaturations}
                                                                                               }};
        const std::array<VectorRef, 2> correspondingPartsOfB = {b.head(systemSize), b.tail(systemSize)};

        const Real estimateOfCapitalB = 1; // see eqn. on page 6194 (last row)
        Real maximumVarianceGrowthFactor = 0;
        Real maximumGamma = 0;
        int maximumVarianceGroupIndex = -1;



        if (doTryToCorrectForGrowth) {
            log()->info("Trying to correct for variance growth.");
        } else {
            log()->info("Not trying to correct for variance growth. Make sure it looks linear!");
        }


        if (onlyCorrectPressurePart) {
            log()->info("Only correcting pressure part");
        } else {
            // dSigma(i)/dSw(i) wouldn't be identity anymore
            std::abort();
        }

        for (int groupIndex = 0; groupIndex < int(jordanBlocksInConsideration.size()); ++groupIndex) {
            VectorRef correspondingPartOfB = correspondingPartsOfB[groupIndex];

            for (int colIndex = 0; colIndex < systemSize; ++colIndex) {
                Real gamma = 0;
                for (SparseMatrix* jordanBlock: jordanBlocksInConsideration[groupIndex]) {
                    gamma += jordanBlock->col(colIndex).cwiseAbs().sum();
                    //log()->debug("gamma += {}, gamma = {}", jordanBlock->col(colIndex).cwiseAbs().sum(), gamma);
                }


                // correct to use 1 - dPi/dp' on the diagonal
                if (groupIndex == 0) {
                    gamma -= std::abs(pressuresByPressures.coeff(colIndex, colIndex));
                    gamma += std::abs(1 - pressuresByPressures.coeff(colIndex, colIndex));

                }

                const Real varianceGrowthFactor =
                      std::pow(gamma, 2) + gamma * std::abs(correspondingPartOfB(colIndex)) / estimateOfCapitalB;

                if (varianceGrowthFactor > maximumVarianceGrowthFactor) {
                    //log()->debug("gamma = {}, bpart = {}, vargrowth = {}", gamma, std::abs(correspondingPartOfB(colIndex)), varianceGrowthFactor);
                    maximumVarianceGrowthFactor = varianceGrowthFactor;
                    maximumVarianceGroupIndex = groupIndex;
                }

                if (gamma > maximumGamma) {
                    maximumGamma = gamma;
                }


                if (doTryToCorrectForGrowth) {
                    const Real correctionQuotient = std::max(1.01, gamma);

                    if (!onlyCorrectPressurePart || groupIndex == 0) {
                        for (SparseMatrix* jordanBlock: jordanBlocksInConsideration[groupIndex]) {
                            jordanBlock->col(colIndex) = jordanBlock->col(colIndex) / correctionQuotient;
                        }
                    }

                    correspondingPartOfB(colIndex) /= correctionQuotient;
                    ASSERT(allFinite(b));
                }



            }
        }

        log()->info("Maximum variance growth factor = {} in group {}", maximumVarianceGrowthFactor, maximumVarianceGroupIndex);
        log()->info("Maximum gamma = {}", maximumGamma);
    }

    saturationsByPressures.makeCompressed();
    pressuresByPressures.makeCompressed();
}
