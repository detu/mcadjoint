#pragma once
#include "typedefs.hpp"
#include "cellindex.hpp"

Matrix computeXDerivative(ConstMatrixRef field, const Real meshWidth);
Matrix computeYDerivative(ConstMatrixRef field, const Real meshWidth);
Real getDerivativeAtCellBorder(CellIndex cell,
                               ConstMatrixRef xDerivative, ConstMatrixRef yDerivative,
                               const CellIndex::Direction whichBorder);
CellIndex borderIndexToCenterIndex(CellIndex borderIndex, const CellIndex::Direction whichBorder);
CellIndex centerIndexToBorderIndex(CellIndex centerIndex, const CellIndex::Direction whichBorder);
Matrix computeTotalDarcyVelocitiesX(ConstMatrixRef totalMobilities, Matrix pressureDerivativesX);
Matrix computeTotalDarcyVelocitiesY(ConstMatrixRef totalMobilities, Matrix pressureDerivativesY);
