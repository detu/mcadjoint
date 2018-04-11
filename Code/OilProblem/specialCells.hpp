//
// Created by Stefano Weidmann on 11.04.18.
//

#pragma once
#include "cellindex.hpp"

inline CellIndex findDrillCell(const int numberOfRows, const int numberOfCols) {
    return {numberOfRows-1, 0};
}

inline CellIndex findWellCell(const int numberOfRows, const int numberOfCols) {
    return {0, numberOfCols-1};
}

inline CellIndex findReferenceCell(const int numberOfRows, const int numberOfCols) {
    return {0, 0};
}
