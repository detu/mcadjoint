#pragma once
#include "randomWalkState.hpp"
#include <vector>
#include "bVectorSurrogate.hpp"
#include "cMatrixSurrogate.hpp"

bool transitionState(RandomWalkState& currentState, const BVectorSurrogate& b,
                     const SparseMatrix& pressureResidualsByPressures,
                     const SparseMatrix& pressureResidualsBySaturationsWater,
                     const SparseMatrix& saturationsWaterResidualsByPressure,
                     const SparseMatrix& saturationsWaterResidualsBySaturationsWater, const int numberOfRows,
                     const int numberOfCols, Rng& rng);
void addNewRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters,
                       const int currentTimelevel, const BVectorSurrogate& b, const CMatrixSurrogate& c,
                       std::vector<RandomWalkState>& candidates, Rng& rng);
void logStatisticsAboutRandomWalks(const std::vector<RandomWalkState>& randomWalks);
