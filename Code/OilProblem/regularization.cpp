//
// Created by Stefano Weidmann on 17.05.18.
//

#include "regularization.hpp"
#include "darcyVelocity.hpp"
#include "pressure.hpp"
#include "fixedParameters.hpp"
#include "logging.hpp"
#include "specialCells.hpp"
#include <Eigen/Dense>
#include "sensitivity.hpp"

//#error "Referenzwert (k - kref)^2 pore volume injected = qDeltat/(Vol * Phi)"
//#error "Adjoint"
// Fixer stimestep
// k_rel = k_test * Delta p rel / delta p test
// Symmetrische Lösung
// Test optimierer nahe optimale Lösung
Real computeRegularizationPenalty(ConstMatrixRef logPermeabilities, const Real referenceLogPermeability, const Real regularizationParameter) {
    return (logPermeabilities.array() - referenceLogPermeability).square().sum();
}

Vector deriveRegularizationPenaltyByLogPermeabilities(ConstMatrixRef logPermeabilities, const Real referenceLogPermeability, const Real regularizationParameter) {
    return 2.0 * regularizationParameter *  (Eigen::Map<const Eigen::Array<Real, Eigen::Dynamic, 1>>(logPermeabilities.data(), logPermeabilities.size(), 1) - referenceLogPermeability);
}

Real computeReferenceLogPermeability(const FixedParameters& params) {
    const int numberOfRows = params.initialSaturationsWater.rows();
    const int numberOfCols = params.initialSaturationsWater.cols();

    const Matrix homogeneousPermeabilities = Matrix::Ones(numberOfRows, numberOfCols);
    const Matrix noWaterSaturations = Matrix::Zero(numberOfRows, numberOfCols);

    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, homogeneousPermeabilities, noWaterSaturations);
    const Vector rhs = computeRhsForPressureSystem(params.inflowPerUnitDepthWater(0), numberOfRows, numberOfCols);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const Vector pressuresAsVector = solvePressurePoissonProblem(pressureSystem, rhs);

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);

    const int computedPressureAtDrill = pressuresAsVector(drillCell.linearIndex(numberOfRows));
    const int measuredPressureAtDrill = params.overPressureDrill(0);

    const Real referencePermeability = computedPressureAtDrill / measuredPressureAtDrill * drillCell(homogeneousPermeabilities);
    log()->info("Reference permeability is {}", referencePermeability);

    const Real referenceLogPermeability = std::log(referencePermeability);

    return referenceLogPermeability;

}

void applyRegularizationIfEnabled(const FixedParameters& params, ConstMatrixRef logPermeabilities, SensitivityAndCost& sensitivityAndCost) {
    if (params.regularizationParameter > 0) {
        log()->info("Regularization enabled");

        const Real referenceLogPermeability = computeReferenceLogPermeability(params);
        const Real regularizationParameter = 1;

        log()->info("Using a regularization parameter of {}", regularizationParameter);

        const Vector regularizationPenalty = deriveRegularizationPenaltyByLogPermeabilities(logPermeabilities, referenceLogPermeability, regularizationParameter);
        log()->info("Regularization penalty norm = {}", regularizationPenalty.norm());
        sensitivityAndCost.sensitivity += regularizationPenalty;

        const Real regularizationPenaltyCost = computeRegularizationPenalty(logPermeabilities, referenceLogPermeability, regularizationParameter);
        log()->info("Regularization penalty cost = {}", regularizationPenaltyCost);
        sensitivityAndCost.cost += regularizationPenaltyCost;


    } else {
        log()->info("Regularization disabled");
    }
}