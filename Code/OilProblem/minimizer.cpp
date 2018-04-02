//
// Created by Stefano Weidmann on 02.04.18.
//

#include "oilproblem.hpp"

MinimizationState doAMinimizerStep(const CostFunction& toMinimize, const MinimizationState& oldState, ConstMatrixRef sensitivity) {
    MinimizationState newState(oldState.parameters.map.rows(), oldState.parameters.map.cols());

    newState.stepSize = oldState.stepSize;

    const Vector gradientDirection = sensitivity * oldState.parameters.vec;
    for (;;) {
        newState.parameters.vec = oldState.parameters.vec - newState.stepSize * gradientDirection;
        newState.cost = toMinimize(newState.parameters.map);

        if (newState.cost < oldState.cost) {
            newState.stepSize += std::max(newState.stepSize / 10, 0.001);
            break;
        }

        newState.stepSize /= 2;
    }

    return newState;


}
