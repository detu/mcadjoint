//
// Created by stefanow on 4/27/18.
//

#include "sensitivity.hpp"
#include <stefCommonHeaders/dev.hpp>
#include "specialCells.hpp"
#include "simulationState.hpp"
#include "forward.hpp"
#include "pressure.hpp"
#include "derivativesForAdjoint.hpp"

Matrix computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities) {
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();
    const int numberOfCells = numberOfParameters;

    // TODO
    Vector pressureRhs(Vector::Zero(numberOfCells));
    SimulationState simulationState(numberOfRows, numberOfCols);

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);




    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, simulationState.saturationsWater);


    std::vector<RandomWalkState> randomWalks;
    Rng rng;

    bool breakthroughHappened = false;
    do {



    } while (!breakthroughHappened);


    const Matrix totalMobilitiesBySaturationsWater = computeTotalMobilitiesDerivedBySaturationsWater(permeabilities, simulationState.saturationsWater, params.dynamicViscosityOil, params.dynamicViscosityWater);
    const SparseMatrix pressureResidualsBySaturationsWater = computePressureResidualsDerivedBySaturationWater(simulationState.pressures, totalMobilities, totalMobilitiesBySaturationsWater);



    Matrix sensitivity;
    return sensitivity;
}


