#pragma once
#include "randomWalkState.hpp"
#include <vector>
#include "bVectorSurrogate.hpp"
#include "cMatrixSurrogate.hpp"

bool transitionState(RandomWalkState& currentState, const BVectorSurrogate& b,
                     const SparseMatrix& pressureResidualsByPressures,
                     const SparseMatrix& pressureResidualsBySaturationsWater,
                     const SparseMatrix& saturationsWaterResidualsByPressure,
                     const SparseMatrix& saturationsWaterResidualsBySaturationsWater,
                     const int numberOfRows, const int numberOfCols, Rng& rng);
std::vector<RandomWalkState> initializeRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters, const BVectorSurrogate& b, const CMatrixSurrogate& c);
void logStatisticsAboutRandomWalks(const std::vector<RandomWalkState>& randomWalks);
