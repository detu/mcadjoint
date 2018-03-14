//
// Created by Stefano Weidmann on 13.03.18.
//

#ifndef MCADJOINT_PROBLEMSELECTION_HPP
#define MCADJOINT_PROBLEMSELECTION_HPP

#include <unordered_map>
#include <string>
#include <cstdint>



enum class Problem: std::int64_t {
    MATCH_FINAL_WITH_INITIAL,
    MATCH_DATA_WITH_INITIAL,
    MATCH_DATA_WITH_VISCOSITY
};

Problem getProblemFromDescription(std::string description);

#endif //MCADJOINT_PROBLEMSELECTION_HPP
