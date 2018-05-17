cmake_minimum_required(VERSION 3.9)
project(OilProblem LANGUAGES CXX)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-int-in-bool-context")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

get_filename_component(thirdPartyDir "${CMAKE_CURRENT_SOURCE_DIR}/../thirdParty" REALPATH)

include("${thirdPartyDir}/stefcommonheaders/CMakeLists.txt")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

option(USE_COTIRE "Use cotire" OFF)

if (USE_COTIRE)
    include(cotire)
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")


set(optimizerIncludeDir "${thirdPartyDir}/CppNumericalSolvers/include")
include_directories("${optimizerIncludeDir}" )

set(eigenIncludeDir "${thirdPartyDir}/eigen")
include_directories("${eigenIncludeDir}")

set(simplematioDir "${thirdPartyDir}/simplematio")
set(simplematioBinaryDir "${CMAKE_BINARY_DIR}/simplematio")
add_subdirectory("${simplematioDir}" "${simplematioBinaryDir}")

set(baseLibHeaders utils.hpp logging.hpp)
set(baseLibSources utils.cpp)
set(oilLibHeaders typedefs.hpp cellindex.hpp regularization.hpp vectorToBeMappedAsMatrix.hpp fixedParameters.hpp simulationState.hpp specialCells.hpp randomWalkState.hpp adjointOptions.hpp preconditioningOptions.hpp
        pressure.hpp darcyVelocity.hpp saturation.hpp forward.hpp derivativesForAdjoint.hpp derivativesForAdjointFD.hpp adjoint.hpp sensitivity.hpp dumpToMatFile.hpp preconditioning.hpp)
set(oilLibSources pressure.cpp regularization.cpp darcyVelocity.cpp saturation.cpp forward.cpp derivativesForAdjoint.cpp derivativesForAdjointFD.cpp adjoint.cpp sensitivity.cpp minimizer.cpp dumpToMatFile.cpp preconditioning.cpp)
set(allCommonHeadersAndSources ${oilLibSources} ${oilLibHeaders} ${baseLibSources} ${baseLibHeaders})


set(quarterFiveSpotTraditionalSources quarterFiveSpotTraditional.cpp)
set(quarterFiveSpotSources quarterFiveSpot.cpp)
set(quarterFiveSpotAdjointSources quarterFiveSpotAdjoint.cpp)
set(testSources test.cpp)


add_executable(Test ${testSources})
add_executable(TestVerbose ${testSources})
add_executable(QuarterFiveSpot ${quarterFiveSpotSources})
add_executable(QuarterFiveSpotAdjoint ${quarterFiveSpotAdjointSources})
target_compile_definitions(TestVerbose PUBLIC -DVERBOSE_TESTS)
add_executable(QuarterFiveSpotTraditional ${quarterFiveSpotTraditionalSources})
add_executable(DerivativesForAdjointTest derivativesForAdjointTest.cpp)

add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int -DEIGEN_NO_AUTOMATIC_RESIZING)
if (lowercaseBuildType STREQUAL "debug")
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN -DEIGEN_DONT_VECTORIZE)
endif()
target_compile_definitions(QuarterFiveSpotTraditional PUBLIC -DJUST_COMPUTE_ADJOINT)

string(TOLOWER CMAKE_BUILD_TYPE lowercaseBuildType)
foreach(executable IN ITEMS Test TestVerbose QuarterFiveSpot QuarterFiveSpotAdjoint QuarterFiveSpotTraditional DerivativesForAdjointTest)
    target_sources("${executable}" PUBLIC "${allCommonHeadersAndSources}")
    target_link_libraries("${executable}" stefCommonHeaders simplematio Threads::Threads)
    if (USE_COTIRE)
        cotire("${executable}")
    endif()
endforeach()

