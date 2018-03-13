//
// Created by Stefano Weidmann on 13.03.18.
//

#ifndef MCADJOINT_PROBLEMSELECTION_HPP
#define MCADJOINT_PROBLEMSELECTION_HPP

#include <map>
#include <string>

enum class Problem {
    MATCH_FINAL_WITH_INITIAL,
    MATCH_DATA_WITH_INITIAL
};

extern const std::map<std::string, Problem> PROBLEM_DESCRIPTIONS;

Problem getProblemFromDescription(const std::string& description);

#endif //MCADJOINT_PROBLEMSELECTION_HPP
