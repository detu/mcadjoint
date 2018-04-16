//
// Created by Stefano Weidmann on 10.04.18.
//

#pragma once

#include "vectorToBeMappedAsMatrix.hpp"
#include "typedefs.hpp"


struct SimulationState {
    Matrix saturationsWater;
    VectorToBeMappedAsMatrix pressures;
    Real time;

    inline SimulationState(const int matrixRows, const int matrixCols):
          saturationsWater(matrixRows, matrixCols), pressures(matrixRows, matrixCols), time(0) {}

    inline SimulationState():
          SimulationState(0, 0) {}



};
