//
// Created by Stefano Weidmann on 05.04.18.
//

#include "logging.hpp"

#include <iostream>
using namespace std;
using namespace spdlog;



static shared_ptr<logger> setUpLog() {

    const char* logFileName = "log";
    const char* logLevelEnv = getenv("LOG_LEVEL");
    if (logLevelEnv == nullptr) {
        logLevelEnv = "off";
    }

    const auto level = level::from_str(logLevelEnv);
    set_level(level);

    cerr << "Logging to file \"" << logFileName << "\" and to stderr. Logging level is " << level::to_str(level) << "\n";


    std::vector<sink_ptr> sinks;
    sinks.push_back(make_shared<sinks::ansicolor_stderr_sink_st>());
    sinks.push_back(make_shared<sinks::simple_file_sink_st>(logFileName));


    return create("logger", sinks.begin(), sinks.end());
}

const shared_ptr<logger> LOGGER = setUpLog();