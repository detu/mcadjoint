cmake_minimum_required(VERSION 3.8)
project(stefCommonHeaders LANGUAGES C CXX VERSION 0.1.0)
option(USE_EIGEN_INITIALIZER_LISTS "Whether to include the eigen initializer list plugin" ON)
if (CMAKE_CURRENT_LIST_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(doBuildTestsByDefault ON)
    message("Building tests for stefcommonheaders")
else()
    message("Not building tests for stefcommonheaders")
    set(doBuildTestsByDefault OFF)
endif()

option(BUILD_STEFCOMMONHEADERS_TESTS "Wheter to build the tests" ${doBuildTestsByDefault})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/thirdParty/cotire/CMake")

function(addEigenAndPluginToTarget target required found_var)
    set(eigenInitializerListDir "${CMAKE_CURRENT_LIST_DIR}/thirdParty/eigenInitializerLists/src")

    if(required OR USE_EIGEN_INITIALIZER_LISTS)
        find_package(Eigen3 REQUIRED)
    else()
        find_package(Eigen3)
    endif()

    set("${found_var}" "${EIGEN3_FOUND}" PARENT_SCOPE)
    if(NOT EIGEN3_FOUND)
        return()
    endif()

    target_include_directories(${target} INTERFACE "${EIGEN3_INCLUDE_DIR}")
    target_include_directories(${target} INTERFACE "${eigenInitializerListDir}")
    target_compile_definitions(${target} INTERFACE
        -DEIGEN_MATRIX_PLUGIN="${eigenInitializerListDir}/Matrix_initializer_list.h"
        -DEIGEN_ARRAY_PLUGIN="${eigenInitializerListDir}/Array_initializer_list.h"
    )

endfunction()

# Wrap it in a function as not to pollute the including CMake with my variables
function(addStefCommonHeaders)
    execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}")

    set(headersDir "${CMAKE_CURRENT_LIST_DIR}/include/stefCommonHeaders")
    set(acutestDir "${CMAKE_CURRENT_LIST_DIR}/thirdParty/acutest/include")
    set(arghDir "${CMAKE_CURRENT_LIST_DIR}/thirdParty/argh")
    set(spdlogDir "${CMAKE_CURRENT_LIST_DIR}/thirdParty/spdlog/include")
    get_filename_component(parentDirOfHeadersDir "${headersDir}" DIRECTORY)


    add_library(stefCommonHeaders INTERFACE)
    target_compile_features(stefCommonHeaders INTERFACE cxx_std_11)
    target_include_directories(stefCommonHeaders INTERFACE "${parentDirOfHeadersDir}" "${acutestDir}" "${arghDir}" "${spdlogDir}")

    addEigenAndPluginToTarget(stefCommonHeaders OFF EIGEN3_FOUND)

    install(TARGETS stefCommonHeaders EXPORT stefCommonHeadersExport)
    export(EXPORT stefCommonHeadersExport FILE "${CMAKE_BINARY_DIR}/stefCommonHeadersConfig.cmake")
    export(PACKAGE stefCommonHeaders)




    if(BUILD_STEFCOMMONHEADERS_TESTS)

        if(EIGEN3_FOUND)
            add_executable(stefCommonHeadersTestBinary test/test.cpp)
            addEigenAndPluginToTarget(stefCommonHeadersTestBinary ON dummyVariableNeverUsed)
            target_link_libraries(stefCommonHeadersTestBinary stefCommonHeaders)
            add_test(NAME stefCommonHeadersTest COMMAND stefCommonHeadersTestBinary)
        endif()

        add_executable(stefFenvTestBinary "${CMAKE_CURRENT_LIST_DIR}/test/stefFenvTest.c")
        target_link_libraries(stefFenvTestBinary stefCommonHeaders m)


        add_executable(enterPositiveNumbers "${CMAKE_CURRENT_LIST_DIR}/test/enterPositiveNumbers.cpp")
        target_link_libraries(enterPositiveNumbers stefCommonHeaders)

        add_executable(testAssert "${CMAKE_CURRENT_LIST_DIR}/test/testAssert.cpp")
        target_link_libraries(testAssert stefCommonHeaders)

        add_test(NAME stefFenvTest COMMAND stefFenvTestBinary)
    endif()
endfunction()
addStefCommonHeaders()
