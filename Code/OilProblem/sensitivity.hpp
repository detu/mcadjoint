#pragma once
#include "typedefs.hpp"
#include "fixedParameters.hpp"

struct SensitivityAndCost {
    Vector sensitivity;
    Real cost;
};

SensitivityAndCost computeSensitivityAndCost(const FixedParameters& params, ConstMatrixRef permeabilities,
                                             ConstMatrixRef logPermeabilities, Rng& rng);

SensitivityAndCost computeSensitivityAndCostTraditional(const FixedParameters& params, ConstMatrixRef permeabilities,
                                                        ConstMatrixRef logPermeabilities);
