//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once
#include "typedefs.hpp"
#include "cellindex.hpp"
#include "pressure.hpp"
#include "utils.hpp"

#ifdef JUST_COMPUTE_ADJOINT
    #pragma message("Computing adjoint in mcmc")
#endif

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


    inline Real operator() (const CellIndex stateCell, const CellIndex parameterCell, const bool isAPressure) const {
        CellIndex derivativeCell = pressureToTransmissibilityIndex(stateCell, parameterCell, numberOfRows);


        #ifdef JUST_COMPUTE_ADJOINT
        if (!isAPressure) {
            derivativeCell.j += endIndexInPressurePart;
        }
        return -Real(derivativeCell.i == derivativeCell.j);
        #else

        if (isAPressure) {
            return derivativeCell(pressureResidualsByLogPermeability);
        } else {
            return derivativeCell(saturationWaterResidualsByLogPermeability);
        }
        #endif
    }

};