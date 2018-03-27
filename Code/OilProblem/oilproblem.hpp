//
// Created by Stefano Weidmann on 24.03.18.
//

#ifndef STEFCOMMONHEADERS_OILPROBLEM_HPP
#define STEFCOMMONHEADERS_OILPROBLEM_HPP

#include "typedefs.hpp"
#include "cellindex.hpp"

#if !defined(__GNUC__) && !defined(__attribute__)
    #define __attribute__(ignored)
#endif

__attribute__((pure))
Real computeTransmissibility(ConstMatrixRef lambdas, const CellIndex& fromCell, const CellIndex& toCell);




__attribute__((pure))
SparseMatrix assembleTransmissibilityMatrix(ConstMatrixRef lambdas);


__attribute__((pure))
VectorToBeMappedAsMatrix solvePressurePoissonProblem(const SparseMatrix& transmissibilities, ConstMatrixRef sources);



#endif //STEFCOMMONHEADERS_OILPROBLEM_HPP
