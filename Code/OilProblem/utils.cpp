//
// Created by Stefano Weidmann on 13.03.18.
//

#include <algorithm>

#include "utils.hpp"

std::string lowerCase(std::string str) {
    transform(str.begin(), str.end(), str.begin(), [] (char c) -> char {return char(tolower(c));});
    return str;
}