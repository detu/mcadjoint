cmake_minimum_required(VERSION 3.8)
project(OilProblem LANGUAGES CXX)

option(USE_MKL_IF_PRESENT "Use Intel MKL if present" OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SIMPLEMATIO REQUIRED simplematio)
link_directories(${SIMPLEMATIO_LIBRARY_DIRS})


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(COTIRE_UNITY_SOURCE_MAXIMUM_NUMBER_OF_INCLUDES -j)
get_filename_component(thirdPartyDir "${CMAKE_CURRENT_SOURCE_DIR}/../thirdParty" REALPATH)
set(COTIRE_ADDITIONAL_PREFIX_HEADER_IGNORE_PATH "${thirdPartyDir}")


include("../thirdParty/stefcommonheaders/CMakeLists.txt")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

set(baseLibHeaders control.h gnuplot_i.h utils.hpp logging.hpp)
set(baseLibSources control.cpp gnuplot_i.cpp utils.cpp logging.cpp)
set(oilLibHeaders oilProblem.hpp typedefs.hpp cellindex.hpp vectorToBeMappedAsMatrix.hpp fixedParameters.hpp simulationState.hpp specialCells.hpp)
set(oilLibSources pressure.cpp darcyVelocity.cpp saturation.cpp minimizer.cpp forward.cpp)
add_library(commonLib STATIC ${oilLibSources} ${oilLibHeaders} ${baseLibSources} ${baseLibHeaders})
target_link_libraries(commonLib stefCommonHeaders)

set(quarterFiveSpotSources quarterFiveSpot.cpp)
set(dryRunSources dryRun.cpp)
set(mainSources main.cpp)
set(testSources test.cpp)


add_executable(DryRun ${dryRunSources})
add_executable(Test ${testSources})
add_executable(TestVerbose ${testSources})
add_executable(QuarterFiveSpot ${quarterFiveSpotSources})
target_compile_definitions(TestVerbose PUBLIC -DVERBOSE_TESTS)
add_executable(OilProblem ${mainSources})

target_link_libraries(QuarterFiveSpot ${SIMPLEMATIO_LIBRARIES})
target_include_directories(QuarterFiveSpot PUBLIC ${SIMPLEMATIO_INCLUDE_DIRS})
target_compile_options(QuarterFiveSpot PUBLIC ${SIMPLEMATIO_CFLAGS_OTHER})

add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int -DEIGEN_NO_AUTOMATIC_RESIZING)
if (lowercaseBuildType STREQUAL "debug")
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
endif()

string(TOLOWER CMAKE_BUILD_TYPE lowercaseBuildType)
foreach(executable IN ITEMS DryRun Test TestVerbose OilProblem QuarterFiveSpot)
    target_link_libraries(${executable} commonLib)
    #cotire(${executable})
endforeach()


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
