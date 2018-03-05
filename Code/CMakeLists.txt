cmake_minimum_required(VERSION 3.8)
project(MCAdjoint LANGUAGES CXX)

add_subdirectory(thirdParty/stefcommonheaders)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

if(NOT (Eigen3_VERSION_MAJOR GREATER 3 OR Eigen3_VERSION_MAJOR EQUAL 3 AND NOT (Eigen3_VERSION_MINOR LESS 3)))
    message(FATAL_ERROR "This package needs Eigen >= 3.3, but found was ${Eigen3_VERSION}!")
endif()
include_directories(. "${EIGEN3_INCLUDE_DIR}")


add_executable(MCAdjoint
        control.cpp
        control.h
        driver.h
        driverBurger.cpp
        gnuplot_i.cpp
        gnuplot_i.h
        main.cpp
        randomgeneral.h
        gnuPlotSetTerminal.cpp
        gnuPlotSetTerminal.hpp)
target_link_libraries(MCAdjoint stefCommonHeaders)

