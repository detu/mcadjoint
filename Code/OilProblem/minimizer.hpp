//
// Created by Stefano Weidmann on 30.04.18.
//

#pragma once

#include "fixedParameters.hpp"

struct PermeabilitiesAndCost {
    Matrix permeabilities;
    Real cost;
};

PermeabilitiesAndCost
matchWithPermeabilities(const FixedParameters& params, const Real tolerance, const int maxIterations);