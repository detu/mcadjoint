//
// Created by Stefano Weidmann on 27.03.18.
//

#ifndef STEFCOMMONHEADERS_PRESSUREINTERNAL_HPP
#define STEFCOMMONHEADERS_PRESSUREINTERNAL_HPP

#include "oilproblem.hpp"

__attribute__((pure))
static inline CellIndex pressureToTransmissibilityIndex(
      const CellIndex& fromCell,
      const CellIndex& toCell,
      const int numberOfRows);

#endif //STEFCOMMONHEADERS_PRESSUREINTERNAL_HPP
