//
// Created by Stefano Weidmann on 13.03.18.
//

#include "problemSelection.hpp"

#include "utils.hpp"
#include <stdexcept>

const std::unordered_map<std::string, Problem> SPARSE_DATA_FROM_DESCRIPTION = {
      {"matchdatawithinitial", Problem::MATCH_DATA_WITH_INITIAL},
      {"matchfinalwithinitial", Problem::MATCH_FINAL_WITH_INITIAL},
      {"matchdatawithviscosity", Problem::MATCH_DATA_WITH_VISCOSITY}
};

Problem getProblemFromDescription(std::string description) {
    description = lowerCase(description);
    try {
        return SPARSE_DATA_FROM_DESCRIPTION.at(description);
    } catch (const std::out_of_range& outOfRange) {
        throw std::invalid_argument("The problem \"" + description + "\" is unknown!");
    }
}