cmake_minimum_required(VERSION 3.10)
project(MCAdjoint)

add_subdirectory(thirdParty/stefcommonheaders)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
include_directories(. "${EIGEN3_INCLUDE_DIR}")


add_executable(MCAdjoint
        control.cpp
        control.h
        driver.h
        driverBurger.cpp
        gnuplot_i.cpp
        gnuplot_i.h
        main.cpp
        randomgeneral.h)
target_link_libraries(MCAdjoint INTERFACE stefCommonHeaders)
