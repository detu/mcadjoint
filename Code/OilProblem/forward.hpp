#pragma once

#include "typedefs.hpp"
#include "fixedParameters.hpp"
#include "simulationState.hpp"
#include "randomWalkState.hpp"

bool stepForwardProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                        SimulationState& currentState);
bool stepForwardAndAdjointProblemTraditional(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                             const int currentTimelevel, SimulationState& simulationState, Matrix& adjointMatrix, Vector& adjointRhs);
bool stepForwardAndAdjointProblem(const FixedParameters& params, ConstMatrixRef permeabilities, const int currentTimelevel,
                                  SimulationState& simulationState, std::vector<RandomWalkState>& randomWalks, std::vector<Rng>& rngs);