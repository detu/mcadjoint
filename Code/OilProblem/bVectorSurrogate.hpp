//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once

#include "typedefs.hpp"
#include "cellindex.hpp"
#include "specialCells.hpp"
#include "derivativesForAdjoint.hpp"
#include "logging.hpp"


struct BVectorSurrogate {
    BVectorSurrogate(const Real computedPressureAtDrill, const Real measuredPressureAtDrill,
                     const int numberOfRows, const int numberOfCols):
          computedPressureAtDrill(computedPressureAtDrill), measuredPressureAtDrill(measuredPressureAtDrill),
          numberOfRows(numberOfRows), endOfPressurePart(numberOfCols * numberOfRows), drillCell(findDrillCell(numberOfRows, numberOfCols)),
          numberOfCols(numberOfCols)
    {
    }

    inline Real operator()(int linearStateIndex) const {
        const bool isAPressure = linearStateIndex < endOfPressurePart;
        if (!isAPressure) {
            linearStateIndex -= endOfPressurePart;
        }
        const CellIndex stateCell = CellIndex::fromLinearIndex(linearStateIndex, numberOfRows);

        return operator()(stateCell, isAPressure);

    }


    inline Real operator()(const CellIndex stateCell, const bool isAPressure) const {
        if (isAPressure) {
            return computeCostFunctionDerivedByPressureEntry(computedPressureAtDrill, measuredPressureAtDrill, numberOfRows, numberOfCols, stateCell);
        } else {
            return computeCostFunctionDerivedBySaturationsWaterEntry(numberOfRows, numberOfCols, stateCell);
        }

    }

    const Real computedPressureAtDrill;
    const Real measuredPressureAtDrill;

    const int endOfPressurePart;
    const int numberOfRows;
    const int numberOfCols;


    const CellIndex drillCell;
};
#pragma GCC poison BVectorSurrogate