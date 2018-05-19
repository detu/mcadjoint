//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once
#include "cellindex.hpp"
struct RandomWalkState {
    CellIndex cell;
    bool isAPressure;
    int currentTimelevel;
    int startingTimeLevel;
    WiderReal W;
    WiderReal D;
    int parameterIndex;
};