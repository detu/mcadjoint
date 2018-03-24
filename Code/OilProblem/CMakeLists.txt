cmake_minimum_required(VERSION 3.10)
project(OilProblem)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include("thirdParty/stefcommonheaders/CMakeLists.txt")
add_executable(OilProblem
        control.cpp
        control.h
        gnuplot_i.cpp
        gnuplot_i.h
        utils.cpp
        utils.hpp oilproblem.cpp oilproblem.hpp typedefs.hpp cellindex.hpp)
target_link_libraries(OilProblem stefCommonHeaders)
