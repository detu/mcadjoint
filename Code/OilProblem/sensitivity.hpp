#pragma once
#include "typedefs.hpp"
#include "fixedParameters.hpp"

Matrix computeSensitivity(const FixedParameters& params, ConstMatrixRef permeabilities);
