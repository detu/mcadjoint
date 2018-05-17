//
// Created by Stefano Weidmann on 17.05.18.
//

#include "pickAnIndex.hpp"
#include <numeric>

int pickAnIndex(const std::vector<Real>& weights, const Real standardUniformRV) {
    const Real weightNormalization = std::accumulate(weights.cbegin(), weights.cend(), 0.0);

    Real cumulativeWeightSum = 0;

    const Real decisionVariable = standardUniformRV * weightNormalization;

    int index = 0;
    for (const Real weight: weights) {
        cumulativeWeightSum += weight;

        if (decisionVariable <= cumulativeWeightSum) {
            return index;
        }
        ++index;
    }

    std::abort();
}