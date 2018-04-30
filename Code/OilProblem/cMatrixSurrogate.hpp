//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once
#include "typedefs.hpp"
#include "cellindex.hpp"
#include "pressure.hpp"

struct CMatrixSurrogate {

    inline CMatrixSurrogate(const SparseMatrix& pressureResidualsByLogPermeability, const SparseMatrix& saturationWaterResidualsByLogPermeability, const int numberOfRows, const int numberOfCols):
          pressureResidualsByLogPermeability(pressureResidualsByLogPermeability),
          saturationWaterResidualsByLogPermeability(saturationWaterResidualsByLogPermeability),
          endIndexInPressurePart(numberOfRows * numberOfCols), numberOfRows(numberOfRows){
    }

    const SparseMatrix& saturationWaterResidualsByLogPermeability;
    const SparseMatrix& pressureResidualsByLogPermeability;
    const int endIndexInPressurePart;
    const int numberOfRows;

    inline Real operator()(int linearStateIndex, const int linearParameterIndex) const {

        const bool isAPressure = linearStateIndex < endIndexInPressurePart;

        if (!isAPressure) {
            linearStateIndex -= endIndexInPressurePart;
        }

        const CellIndex stateCell = CellIndex::fromLinearIndex(linearStateIndex, numberOfRows);
        const CellIndex parameterCell = CellIndex::fromLinearIndex(linearParameterIndex, numberOfRows);


    }

    inline Real operator() (const CellIndex stateCell, const CellIndex parameterCell, const bool isAPressure) const {
        const CellIndex derivativeCell = pressureToTransmissibilityIndex(stateCell, parameterCell, numberOfRows);
        if (isAPressure) {
            return derivativeCell(pressureResidualsByLogPermeability);
        } else {
            return derivativeCell(saturationWaterResidualsByLogPermeability);
        }
    }
};