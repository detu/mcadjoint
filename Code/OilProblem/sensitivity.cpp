//
// Created by stefanow on 4/27/18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/dev.hpp>

#ifdef TODO
Matrix computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities) {
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();
    const int numberOfCells = numberOfParameters;

    // TODO
    Vector press

    SimulationState simulationState(numberOfRows, numberOfCols);
    stepForwardProblem(params, permeabilities, simulationState)
    const BVectorSurrogate b()

    initializeRandomWalks(numberOfRows, numberOfCols, numberOfParameters,


    return sensitivity;
}

Real computeCost(const FixedParameters& params, ConstMatrixRef permeabilities) {
    const int numberOfPermeabilities = permea
    initializeRandomWalks(params.numberOfRows, params.numberOfCols)
    // TODO

    DEV_STUB();
    return NAN;


}


CostFunction getCostFunctionForMinimizer(const FixedParameters& params) {
    return [&] (ConstMatrixRef logPermeabilities) {
        return computeCost(params, logPermeabilities.array().exp().matrix());
    };
}
#endif
