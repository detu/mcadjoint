//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once

#include "typedefs.hpp"
#include "cellindex.hpp"



struct BVectorSurrogate {
    BVectorSurrogate(const Real computedPressureAtDrill, const Real measuredPressureAtDrill,
                     const int numberOfRows, const int numberOfCols):
          computedPressureAtDrill(computedPressureAtDrill), measuredPressureAtDrill(measuredPressureAtDrill),
          numberOfRows(numberOfRows), endOfPressurePart(numberOfCols * numberOfRows), drillCell(findDrillCell(numberOfRows, numberOfCols)) {}

    inline Real operator()(const int linearStateIndex) const {
        const bool isAPressure = linearStateIndex >= endOfPressurePart;

        if (!isAPressure) {
            return 0;
        }

        const CellIndex stateCell = CellIndex::fromLinearIndex(linearStateIndex, numberOfRows);

        return operator()(stateCell, isAPressure);
    }


    inline Real operator()(const CellIndex stateCell, const bool isAPressure) const {
        if (!isAPressure) {
            return 0;
        }


        if (stateCell == drillCell) {
            return 2 * (computedPressureAtDrill - measuredPressureAtDrill);
        } else {
            return 0;
        }
    }
    const Real computedPressureAtDrill;
    const Real measuredPressureAtDrill;

    const int endOfPressurePart;
    const int numberOfRows;

    const CellIndex drillCell;
};