cmake_minimum_required(VERSION 3.10)
project(OilProblem LANGUAGES CXX)

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

option(USE_MKL_IF_PRESENT "Use Intel MKL if present" OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SIMPLEMATIO REQUIRED simplematio)
link_directories(${SIMPLEMATIO_LIBRARY_DIRS})


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

get_filename_component(thirdPartyDir "${CMAKE_CURRENT_SOURCE_DIR}/../thirdParty" REALPATH)

include("../thirdParty/stefcommonheaders/CMakeLists.txt")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

set(useCotire OFF)

if (useCotire)
    set(COTIRE_UNITY_SOURCE_MAXIMUM_NUMBER_OF_INCLUDES -j)
    set(COTIRE_ADDITIONAL_PREFIX_HEADER_IGNORE_PATH "${thirdPartyDir}")
    set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY COTIRE_ENABLE_PRECOMPILED_HEADER FALSE)
    include(cotire)
endif()


set(baseLibHeaders control.h gnuplot_i.h utils.hpp logging.hpp)
set(baseLibSources control.cpp gnuplot_i.cpp utils.cpp logging.cpp)
set(oilLibHeaders oilProblem.hpp typedefs.hpp cellindex.hpp vectorToBeMappedAsMatrix.hpp fixedParameters.hpp simulationState.hpp specialCells.hpp)
set(oilLibSources pressure.cpp darcyVelocity.cpp saturation.cpp minimizer.cpp forward.cpp derivativesForAdjoint.cpp adjoint.cpp adjointState.hpp bVectorSurrogate.hpp cMatrixSurrogate.hpp randomWalkState.hpp)
add_library(commonLib STATIC ${oilLibSources} ${oilLibHeaders} ${baseLibSources} ${baseLibHeaders})
target_link_libraries(commonLib stefCommonHeaders)


set(quarterFiveSpotSources quarterFiveSpot.cpp)
set(quarterFiveSpotUnitySources ${quarterFiveSpotSources} ${oilLibSources} ${oilLibHeaders} ${baseLibSources} ${baseLibHeaders})
set(dryRunSources dryRun.cpp)
set(mainSources main.cpp)
set(testSources test.cpp)


add_executable(DryRun ${dryRunSources})
add_executable(Test ${testSources})
add_executable(TestVerbose ${testSources})
add_executable(QuarterFiveSpot ${quarterFiveSpotSources})
target_compile_definitions(TestVerbose PUBLIC -DVERBOSE_TESTS)
add_executable(OilProblem ${mainSources})


foreach(target IN ITEMS QuarterFiveSpot)
    target_link_libraries(${target} ${SIMPLEMATIO_LIBRARIES})
    target_include_directories(${target} PUBLIC ${SIMPLEMATIO_INCLUDE_DIRS})
    target_compile_options(${target} PUBLIC ${SIMPLEMATIO_CFLAGS_OTHER})
endforeach()

add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int -DEIGEN_NO_AUTOMATIC_RESIZING)
if (lowercaseBuildType STREQUAL "debug")
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
endif()

string(TOLOWER CMAKE_BUILD_TYPE lowercaseBuildType)
foreach(executable IN ITEMS DryRun Test TestVerbose OilProblem QuarterFiveSpot)
    target_link_libraries(${executable} commonLib)
    if (useCotire)
        cotire(${executable})
    endif()
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
