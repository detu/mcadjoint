//
// Created by Stefano Weidmann on 02.04.18.
//

#pragma once

#include "typedefs.hpp"

struct FixedParameters {
    Real dynamicViscosityOil;
    Real dynamicViscosityWater;
    Real porosity;

    Real finalTime;

    Real meshWidth;

    int maxNumberOfTimesteps;

    Real regularizationParameter;


    WellFunction inflowPerUnitDepthWater;

    PressureFunction overPressureDrill;

    Matrix initialSaturationsWater;
    Matrix initialPermeabilities;

    int numberOfRandomWalksToAdd;

    bool enableAntitheticSampling;
    bool symmetrizeGradient;


};
