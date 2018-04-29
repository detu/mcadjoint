#pragma once
#include "typedefs.hpp"
#include "fixedParameters.hpp"

Vector computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities);
