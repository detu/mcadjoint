//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once
#include "cellindex.hpp"
struct RandomWalkState {
    CellIndex cell;
    bool isAPressure;
    int currentTimelevel;
    long double W;
    long double D;
    int parameterIndex;
};