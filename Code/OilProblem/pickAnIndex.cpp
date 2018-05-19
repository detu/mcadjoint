//
// Created by Stefano Weidmann on 17.05.18.
//

#include "pickAnIndex.hpp"
#include <numeric>
#include <algorithm>

constexpr bool useQuantileBasedApproach = true;

int pickAnIndex(const std::vector<Real>& weights, const Real standardUniformRV) {
    if (useQuantileBasedApproach) {
        std::vector<Real> weightsForQuantile = weights;
        auto nthIterator = weightsForQuantile.begin() + int(standardUniformRV * weightsForQuantile.size());
        if (nthIterator == weightsForQuantile.end()) {
            --nthIterator;
        }
        std::nth_element(weightsForQuantile.begin(), nthIterator, weightsForQuantile.end());

        const int indexInOriginalWeights = std::find(weights.cbegin(), weights.cend(), *nthIterator) - weights.cbegin();
        return indexInOriginalWeights;
    } else {
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


}