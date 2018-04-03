//
// Created by Stefano Weidmann on 02.04.18.
//

#ifndef STEFCOMMONHEADERS_FIXEDPARAMETERS_HPP
#define STEFCOMMONHEADERS_FIXEDPARAMETERS_HPP

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

#endif //STEFCOMMONHEADERS_FIXEDPARAMETERS_HPP
