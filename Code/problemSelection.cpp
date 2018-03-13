//
// Created by Stefano Weidmann on 13.03.18.
//

#include "problemSelection.hpp"

#include <algorithm>
#include <stdexcept>

const std::map<std::string, Problem> PROBLEM_DESCRIPTIONS = {
      {"matchdatawithinitial", Problem::MATCH_DATA_WITH_INITIAL},
      {"matchfinalwithinitial", Problem::MATCH_FINAL_WITH_INITIAL}
};

Problem getProblemFromDescription(std::string description) {
    std::transform(description.begin(), description.end(), description.begin(), std::tolower);
    try {
        return PROBLEM_DESCRIPTIONS.at(description);
    } catch (const std::out_of_range& outOfRange) {
        throw std::invalid_argument("The problem \"" + description + "\" is unknown!");
    }
}