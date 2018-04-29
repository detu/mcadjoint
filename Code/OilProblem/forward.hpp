#pragma once

#include "typedefs.hpp"
#include "fixedParameters.hpp"
#include "simulationState.hpp"
#include "randomWalkState.hpp"

bool stepForwardProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                        SimulationState& currentState);
void stepForwardAndAdjointProblem(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& simulationState, std::vector<RandomWalkState>& randomWalks, Rng& rng);