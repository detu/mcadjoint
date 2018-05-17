#pragma once
#include "randomWalkState.hpp"
#include <list>

bool transitionState(RandomWalkState& currentState, ConstVectorRef b,
                     const SparseMatrix& pressureResidualsByPressures,
                     const SparseMatrix& pressureResidualsBySaturationsWater,
                     const SparseMatrix& saturationsWaterResidualsByPressure,
                     const SparseMatrix& saturationsWaterResidualsBySaturationsWater, const int numberOfRows,
                     const int numberOfCols, const Real standardUniformNumber);

void addNewRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters,
                       const int currentTimelevel, ConstVectorRef b, SparseMatrix c,
                       std::list<RandomWalkState>& randomWalks, std::list<RandomWalkState>& antitheticRandomWalks, Rng& rng);

void logStatisticsAboutRandomWalks(const std::list<RandomWalkState>& randomWalks);



void removeAbsorbedStates(std::list<RandomWalkState>& randomWalks,
                          Eigen::VectorXi& numberOfRemovedAbsorbedStates, Vector& sumOfDValuesOfAbsorbedStates);