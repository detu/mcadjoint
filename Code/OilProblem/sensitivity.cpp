//
// Created by stefanow on 4/27/18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/dev.hpp>




Matrix computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities) {
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();
    const int numberOfCells = numberOfParameters;

    // TODO
    Vector pressureRhs(Vector::Zero(numberOfCells));
    SimulationState simulationState(numberOfRows, numberOfCols);

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);

    // initialization
    std::vector<RandomWalkState> randomWalks = initializeRandomWalks(params, permeabilities, simulationState);
    std::vector<RandomWalkState> stuckRandomWalks;

    do {

        do {

        } while (!stuckRandomWalks.empty());
    } while (!breakthroughHappened);


    const Matrix totalMobilitiesBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    const SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures, totalMobilities, totalMobilitiesBySaturationsWater);




    return sensitivity;
}



CostFunction getCostFunctionForMinimizer(const FixedParameters& params) {
    return [&] (ConstMatrixRef logPermeabilities) {
        return computeCost(params, logPermeabilities.array().exp().matrix());
    };
}