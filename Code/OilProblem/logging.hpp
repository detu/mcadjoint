//
// Created by Stefano Weidmann on 05.04.18.
//

#pragma once
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <memory>


#ifdef MULTITHREADING
#include <stefCommonHeaders/omp_mutex.hpp>
extern omp_mutex LOG_MUTEX;
#endif


inline std::shared_ptr<spdlog::logger> log() {
    #ifdef MULTITHREADING
    std::lock_guard<omp_mutex> lockGuard(LOG_MUTEX);
    #endif

    return spdlog::get("logger");

}

