//
// Created by Stefano Weidmann on 17.05.18.
//

#pragma once
#include "typedefs.hpp"
#include "fixedParameters.hpp"
#include "sensitivity.hpp"

Real computeRegularizationPenalty(ConstMatrixRef logPermeabilities, const Real referenceLogPermeability, const Real regularizationParameter);
Real computeReferenceLogPermeability(const FixedParameters& params);
Vector deriveRegularizationPenaltyByLogPermeabilities(ConstMatrixRef logPermeabilities, const Real referenceLogPermeability, const Real regularizationParameter);
void applyRegularizationIfEnabled(const FixedParameters& params, ConstMatrixRef logPermeabilities, SensitivityAndCost& sensitivityAndCost);
