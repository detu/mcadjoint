//
// Created by Stefano Weidmann on 05.04.18.
//

#include "logging.hpp"

#include <stefCommonHeaders/logging.hpp>


using namespace std;
using namespace spdlog;

const shared_ptr<logger> LOGGER = ::stefCommonHeaders::setUpLog();