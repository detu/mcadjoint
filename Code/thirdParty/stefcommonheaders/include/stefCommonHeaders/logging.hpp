//
// Created by Stefano Weidmann on 08.04.18.
//

#ifndef STEFCOMMONHEADERS_LOGGING_HPP
#define STEFCOMMONHEADERS_LOGGING_HPP


#include <cstdlib>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/null_sink.h>
#include <iostream>
#include <string>

namespace stefCommonHeaders {

    inline spdlog::level::level_enum getLevelFromEnvVar(const char* envVarName = "LOG_LEVEL") {
        using namespace spdlog;
        const char* logLevelEnv = getenv(envVarName);
        if (logLevelEnv == nullptr) {
            logLevelEnv = "debug";
        }


        const spdlog::level::level_enum level = level::from_str(logLevelEnv);

        return level;
    }

    inline std::shared_ptr<spdlog::logger> setUpNullLog() {
        return std::make_shared<spdlog::logger>("null logger", std::make_shared<spdlog::sinks::null_sink_st>());
    }

    using NoMutex = spdlog::details::null_mutex;

    template <typename Mutex>
    std::shared_ptr<spdlog::logger> setUpLog(const spdlog::level::level_enum level, const bool doLogToStderr = true, const char* logFileName = "log.log") {
        using namespace std;
        using namespace spdlog;

        constexpr bool truncate = true;


        cerr << "------ LOGGING INFO ------\n";



        const bool doLogToFile = logFileName != nullptr && *logFileName != '\0';

        const bool loggingLevelIsOff = level == level::off;

        std::vector<sink_ptr> sinks;
        if (loggingLevelIsOff || (!doLogToFile && !(doLogToStderr))) {
            cerr << "Not logging at all, because ";

            if (loggingLevelIsOff) {
                cerr << "logging off.\n";
            } else {
                cerr << "neither logging to file (nullptr or empty file name given) nor to stderr.\n";
            }
        } else {

            if (doLogToFile) {
                cerr << "Logging to file \"" << logFileName << "\".\n";
                sinks.push_back(make_shared<sinks::simple_file_sink<Mutex>>(logFileName, truncate));
            } else {
                cerr << "Not logging to a file.\n";
            }

            if (doLogToStderr) {
                cerr << "Logging to stderr.\n";
                sinks.push_back(make_shared<sinks::ansicolor_stderr_sink<Mutex>>());
            } else {
                cerr << "Not logging to stderr.\n";
            }

            cerr << "Logging level is \"" << level::to_str(level) << "\".\n";

        }

        cerr << "--------------------------\n";




        shared_ptr<logger> logger = details::registry_t<Mutex>::instance().create("logger", sinks.begin(), sinks.end());
        logger->set_level(level);

        return logger;
    }
}

#endif //STEFCOMMONHEADERS_LOGGING_HPP
