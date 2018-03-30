cmake_minimum_required(VERSION 3.8)
project(OilProblem LANGUAGES CXX)

option(USE_MKL_IF_PRESENT "Use Intel MKL if present" OFF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include("thirdParty/stefcommonheaders/CMakeLists.txt")
include(cotire)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

set(baseLibHeaders control.h gnuplot_i.h utils.hpp)
set(baseLibSources control.cpp gnuplot_i.cpp utils.cpp)
set(oilLibHeaders oilProblem.hpp typedefs.hpp cellindex.hpp)
set(oilLibSources pressure.cpp darcyVelocity.cpp saturation.cpp)
set(dryRunSources dryRun.cpp)
set(mainSources main.cpp)
set(testSources test.cpp)

add_library(BaseLib STATIC ${baseLibSources} ${baseLibHeaders})
target_link_libraries(BaseLib stefCommonHeaders)

add_library(OilLib STATIC ${oilLibSources} ${oilLibHeaders})
target_link_libraries(OilLib BaseLib)

add_executable(DryRun ${dryRunSources})
add_executable(Test ${testSources})
add_executable(OilProblem ${mainSources})
target_link_libraries(DryRun OilLib)
target_link_libraries(OilProblem OilLib)
target_link_libraries(Test OilLib)


if (USE_MKL_IF_PRESENT)
    message("Using MKL enabled by user")
    set(BUILD_SHARED_LIBS OFF)
    find_package(MKL)
    if (MKL_FOUND)
        message("MKL found, using that")
        include_directories(${MKL_INCLUDE_DIRS})
        link_directories(${MKL_LIB_DIRS})
        target_link_libraries(DryRun ${MKL_LIBRARIES})
        target_link_libraries(OilProblem ${MKL_LIBRARIES})
        target_link_libraries(Test ${MKL_LIBRARIES})
        add_definitions(-DEIGEN_USE_MKL_ALL -DUSE_PARDISO)
    else()
        message("MKL not found, using builtin eigen")
    endif()
else()
    message("Using MKL disabled")
endif()
