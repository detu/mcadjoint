cmake_minimum_required(VERSION 3.9)
project(OilProblem LANGUAGES CXX)

set(useOpenMP OFF)

if (useOpenMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    add_definitions(-DMULTITHREADED)
endif()

option(USE_MKL_IF_PRESENT "Use Intel MKL if present" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

get_filename_component(thirdPartyDir "${CMAKE_CURRENT_SOURCE_DIR}/../thirdParty" REALPATH)

include("${thirdPartyDir}/stefcommonheaders/CMakeLists.txt")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

set(useCotire ON)

if (useCotire)
    set(COTIRE_UNITY_SOURCE_MAXIMUM_NUMBER_OF_INCLUDES -j)
    set(COTIRE_ADDITIONAL_PREFIX_HEADER_IGNORE_PATH "${thirdPartyDir}")
    #set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY COTIRE_ENABLE_PRECOMPILED_HEADER FALSE)
    include(cotire)
endif()

set(bgfsIncludeDir "${thirdPartyDir}/LBFGSpp/include")
include_directories("${bgfsIncludeDir}")

set(eigenIncludeDir "${thirdPartyDir}/eigen")
include_directories("${eigenIncludeDir}")

set(simplematioDir "${thirdPartyDir}/simplematio")
set(simplematioBinaryDir "${CMAKE_BINARY_DIR}/simplematio")
add_subdirectory("${simplematioDir}" "${simplematioBinaryDir}")

set(baseLibHeaders control.h gnuplot_i.h utils.hpp logging.hpp)
set(baseLibSources control.cpp gnuplot_i.cpp utils.cpp logging.cpp)
set(oilLibHeaders typedefs.hpp cellindex.hpp vectorToBeMappedAsMatrix.hpp fixedParameters.hpp simulationState.hpp specialCells.hpp bVectorSurrogate.hpp cMatrixSurrogate.hpp randomWalkState.hpp
        pressure.hpp darcyVelocity.hpp saturation.hpp forward.hpp derivativesForAdjoint.hpp adjoint.hpp sensitivity.hpp dumpToMatFile.hpp)
set(oilLibSources pressure.cpp darcyVelocity.cpp saturation.cpp forward.cpp derivativesForAdjoint.cpp adjoint.cpp sensitivity.cpp minimizer.cpp dumpToMatFile.cpp)
add_library(commonLib STATIC ${oilLibSources} ${oilLibHeaders} ${baseLibSources} ${baseLibHeaders})
target_link_libraries(commonLib stefCommonHeaders simplematio)

if(useCotire)
  cotire(commonLib)
endif()

set(quarterFiveSpotTraditionalSources quarterFiveSpotTraditional.cpp)
set(quarterFiveSpotSources quarterFiveSpot.cpp)
set(quarterFiveSpotAdjointSources quarterFiveSpotAdjoint.cpp)
set(mainSources main.cpp)
set(testSources test.cpp)


add_executable(Test ${testSources})
add_executable(TestVerbose ${testSources})
add_executable(QuarterFiveSpot ${quarterFiveSpotSources})
add_executable(QuarterFiveSpotAdjoint ${quarterFiveSpotAdjointSources})
target_compile_definitions(TestVerbose PUBLIC -DVERBOSE_TESTS)
add_executable(OilProblem ${mainSources})
add_executable(QuarterFiveSpotTraditional ${quarterFiveSpotTraditionalSources})


add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int -DEIGEN_NO_AUTOMATIC_RESIZING)
if (lowercaseBuildType STREQUAL "debug")
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)
endif()

string(TOLOWER CMAKE_BUILD_TYPE lowercaseBuildType)
foreach(executable IN ITEMS Test TestVerbose OilProblem QuarterFiveSpot QuarterFiveSpotAdjoint QuarterFiveSpotTraditional)
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
        target_link_libraries(OilProblem ${MKL_LIBRARIES})
        target_link_libraries(Test ${MKL_LIBRARIES})
        add_definitions(-DEIGEN_USE_MKL_ALL -DUSE_PARDISO)
    else()
        message("MKL not found, using builtin eigen")
    endif()
else()
    message("Using MKL disabled")
endif()
