//
// Created by Stefano Weidmann on 13.03.18.
//

#ifndef MCADJOINT_SPARSEDATASELECTION_HPP
#define MCADJOINT_SPARSEDATASELECTION_HPP

#include <string>

enum class SparseData {
    ALWAYS_ZERO,
    HAT_PATTERN,
    ZERO_DIAGONAL,
    LINEAR_OFFDIAGONAL
};

SparseData getSparseDataFromDescription(std::string description);

#endif //MCADJOINT_SPARSEDATASELECTION_HPP
