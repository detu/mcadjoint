//
// Created by Stefano Weidmann on 13.05.18.
//

#include "derivativesForAdjointFD.hpp"
#include "specialCells.hpp"
#include "pressure.hpp"
#include <array>







Real computePressureResidualFDEntry(
      MatrixRef pressures, MatrixRef saturations, MatrixRef logPermeabilities, const CellIndex& cell, const Shift& shift, const FixedParameters& params
) {

    Matrix permeabilities = logPermeabilities.array().exp().matrix();

    std::array<MatrixRef, 3> matricesToShift = {pressures, saturations, Eigen::Ref<Matrix>(permeabilities)};

    Real& valueToShift = shift.cell(matricesToShift[int(shift.where)]);

    const int numberOfRows = pressures.rows();
    const int numberOfCols = pressures.cols();

    const CellIndex wellCell = findWellCell(numberOfRows, numberOfCols);

    const auto computePressureResidualEntry = [&wellCell, &cell, numberOfRows, numberOfCols, &params] (ConstMatrixRef pressures, ConstMatrixRef saturations, ConstMatrixRef permeabilities) -> Real {
        if (cell == wellCell) {
            return wellCell(pressures);
        } else {
            Real residual = 0;
            for (const auto& neighbor: cell.neighbors(numberOfRows, numberOfCols)) {
                residual += computeTransmissibility(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, saturations, cell, neighbor) * (cell(pressures) - neighbor(pressures));
            }

            return residual;
        }
    };


    const Real unshiftedValue = valueToShift;

    if (shift.where == Shift::ShiftWhere::LOG_PERMEABILITIES) {
        valueToShift *= std::exp(shift.amount/2);
    } else {
        valueToShift += shift.amount/2;
    }

    const Real upShiftedResidualValue = computePressureResidualEntry(pressures, saturations, permeabilities);
    valueToShift = unshiftedValue;


    if (shift.where == Shift::ShiftWhere::LOG_PERMEABILITIES) {
        valueToShift /= std::exp(shift.amount/2);
    } else {
        valueToShift -= shift.amount/2;
    }

    const Real downshiftedResidualValue = computePressureResidualEntry(pressures, saturations, permeabilities);

    valueToShift = unshiftedValue;




    return (upShiftedResidualValue - downshiftedResidualValue) / shift.amount;
}

Matrix derivePressureResidualsWithFiniteDifferences(Matrix pressures, Matrix saturations, Matrix logPermeabilities, const Shift::ShiftWhere derivedBy, const FixedParameters& params) {
    const int numberOfRows = pressures.rows();
    const int numberOfCols = pressures.cols();
    const int numberOfPairs = numberOfRows * numberOfCols;
    Matrix derivatives(numberOfPairs, numberOfPairs);

    CellIndex cell = {0, 0};
    CellIndex shiftCell = {0, 0};

    for (cell.j = 0; cell.j < numberOfCols; ++cell.j) {
        for (cell.i = 0; cell.i < numberOfRows; ++cell.i) {
            for (shiftCell.j = 0; shiftCell.j < numberOfCols; ++shiftCell.j) {
                for (shiftCell.i = 0; shiftCell.i < numberOfRows; ++shiftCell.i) {
                    const CellIndex meToShift = pressureToTransmissibilityIndex(cell, shiftCell, numberOfRows);
                    const Shift shift = {shiftCell, derivedBy, 1e-8};
                    meToShift(derivatives) = computePressureResidualFDEntry(pressures, saturations, logPermeabilities, cell, shift, params);
                }
            }
        }
    }

    return derivatives;

}

