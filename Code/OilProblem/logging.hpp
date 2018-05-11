//
// Created by Stefano Weidmann on 05.04.18.
//

#pragma once
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <memory>


inline std::shared_ptr<spdlog::logger> log() {
    return spdlog::get("logger");
}

