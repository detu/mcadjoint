//
// Created by Stefano Weidmann on 02.04.18.
//

#pragma once

#include "typedefs.hpp"

struct FixedParameters {
    Real dynamicViscosityOil;
    Real dynamicViscosityWater;

    WellFunction outFlowPerUnitDepthOil;
    WellFunction outFlowPerUnitDepthWater;

    PressureFunction pressureDrill;
    PressureFunction pressureWell;

    Real meshWidth;

};
