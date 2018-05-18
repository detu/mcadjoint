#pragma once

#include "typedefs.hpp"
#include "fixedParameters.hpp"
#include "simulationState.hpp"
#include "randomWalkState.hpp"
#include <list>

bool stepForwardProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                        SimulationState& currentState);
bool stepForwardAndAdjointProblemTraditional(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                             const int currentTimelevel, SimulationState& simulationState, Matrix& adjointMatrix, Vector& adjointRhs, Matrix& completeC);
bool stepForwardAndAdjointProblem(const FixedParameters& params, const Eigen::Ref<const Matrix>& permeabilities,
                                  const int currentTimelevel, const int numberOfRandomWalksToAdd, SimulationState& simulationState,
                                  std::list<RandomWalkState>& randomWalks, std::list<RandomWalkState>& antitheticRandomWalks, Rng& rng);