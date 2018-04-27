//
// Created by stefanow on 4/27/18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/dev.hpp>

Matrix computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities) {
    Matrix sensitivity;



    return sensitivity;
}

Real computeCost(const FixedParameters& params, ConstMatrixRef permeabilities) {
    DEV_STUB();
    // TODO

    return NAN;
}


CostFunction getCostFunctionForMinimizer(const FixedParameters& params) {
    return [&] (ConstMatrixRef logPermeabilities) {
        return computeCost(params, logPermeabilities.array().exp().matrix());
    };
}
