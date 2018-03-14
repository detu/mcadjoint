//
// Created by Stefano Weidmann on 13.03.18.
//

#include "sparseDataSelection.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include "utils.hpp"

const std::unordered_map<std::string, SparseData> SPARSE_DATA_FROM_DESCRIPTION = {
      {"alwayszero", SparseData::ALWAYS_ZERO},
      {"hatpattern", SparseData::HAT_PATTERN},
      {"zerodiagonal", SparseData::ZERO_DIAGONAL}
};

SparseData getSparseDataFromDescription(std::string description) {
    description = lowerCase(description);
    try {
        return SPARSE_DATA_FROM_DESCRIPTION.at(description);
    } catch (const std::out_of_range& outOfRange) {
        throw std::invalid_argument("The sparse data \"" + description + "\" is unknown!");
    }
}