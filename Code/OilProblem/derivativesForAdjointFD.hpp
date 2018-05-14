//
// Created by Stefano Weidmann on 13.05.18.
//

#pragma once
#include "typedefs.hpp"
#include "fixedParameters.hpp"
#include "cellindex.hpp"


struct Shift {
    enum class ShiftWhere {
        PRESSURES = 0,
        SATURATIONS = 1,
        LOG_PERMEABILITIES = 2
    };

    const CellIndex cell;
    const ShiftWhere where;
    const Real amount;
};

Matrix derivePressureResidualsWithFiniteDifferences(Matrix pressures, Matrix saturations, Matrix logPermeabilities, const Shift::ShiftWhere derivedBy, const FixedParameters& params);