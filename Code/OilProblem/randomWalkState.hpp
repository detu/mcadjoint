//
// Created by Stefano Weidmann on 24.04.18.
//

#pragma once
#include "cellindex.hpp"
#include <memory>
struct RandomWalkState {
    CellIndex cell;
    bool isAPressure;
    int currentTimelevel;
    Real W;
    Real D;
    int parameterIndex;
    constexpr static bool enableAntitheticRandomWalks = true;
    std::unique_ptr<RandomWalkState> antitheticRandomWalk;
    #error "TODO antithetic states, put constexpr into a separate header"
};