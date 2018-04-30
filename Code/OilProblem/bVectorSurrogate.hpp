//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once

#include "typedefs.hpp"
#include "cellindex.hpp"
#include "specialCells.hpp"
#include "derivativesForAdjoint.hpp"


struct BVectorSurrogate {
    BVectorSurrogate(const Real computedPressureAtDrill, const Real measuredPressureAtDrill,
                     const int numberOfRows, const int numberOfCols):
          computedPressureAtDrill(computedPressureAtDrill), measuredPressureAtDrill(measuredPressureAtDrill),
          numberOfRows(numberOfRows), endOfPressurePart(numberOfCols * numberOfRows), drillCell(findDrillCell(numberOfRows, numberOfCols)),
         costFunctionByPressure(computeCostFunctionDerivedByPressure(computedPressureAtDrill, measuredPressureAtDrill, numberOfRows, numberOfCols)),
         costFunctionBySaturationsWater(computeCostFunctionDerivedBySaturationsWater(numberOfRows, numberOfCols))
    {}

    inline Real operator()(int linearStateIndex) const {
        const bool isAPressure = linearStateIndex < endOfPressurePart;
        const CellIndex stateCell = CellIndex::fromLinearIndex(linearStateIndex, numberOfRows);

        const SparseVector& correspondingPartOfBVector = isAPressure? costFunctionByPressure: costFunctionBySaturationsWater;

        if (!isAPressure) {
            linearStateIndex -= endOfPressurePart;
        }

        return correspondingPartOfBVector.coeff(linearStateIndex);

    }


    inline Real operator()(const CellIndex stateCell, const bool isAPressure) const {

        const int linearIndex = stateCell.linearIndex(numberOfRows);

        const SparseVector& correspondingPartOfBVector = isAPressure? costFunctionByPressure: costFunctionBySaturationsWater;

        return correspondingPartOfBVector.coeff(linearIndex);
    }

    const Real computedPressureAtDrill;
    const Real measuredPressureAtDrill;

    const int endOfPressurePart;
    const int numberOfRows;
    const SparseVector costFunctionByPressure;
    const SparseVector costFunctionBySaturationsWater;


    const CellIndex drillCell;
};