//
// Created by Stefano Weidmann on 17.05.18.
//

#pragma once
#include "typedefs.hpp"

Real computeRegularizationPenalty(ConstMatrixRef logPermeabilities, const Real meshWidth);

Vector deriveRegularizationPenaltyByLogPermeabilities(ConstMatrixRef logPermeabilities, const Real meshWidth);