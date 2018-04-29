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

Vector computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities) {
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();
    const int numberOfCells = numberOfParameters;

    SimulationState simulationState(numberOfRows, numberOfCols);

    std::vector<RandomWalkState> randomWalks;
    Rng rng;
    bool breakthroughHappened = false;
    do {
        breakthroughHappened = stepForwardAndAdjointProblem(params, permeabilities, simulationState, randomWalks, rng);
    } while (!breakthroughHappened && simulationState.time < params.finalTime);



    Vector sensitivity(Vector::Zero(numberOfParameters));
    Vector numberOfRandomWalksPerParameter(Vector::Zero(numberOfParameters));


    for (const RandomWalkState& randomWalk: randomWalks) {
        sensitivity(randomWalk.parameterIndex) += randomWalk.D;
        ++numberOfRandomWalksPerParameter(randomWalk.parameterIndex);
    }

    sensitivity.array() /= numberOfRandomWalksPerParameter.array();
    return sensitivity;
}


