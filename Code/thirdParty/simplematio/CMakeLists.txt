cmake_minimum_required(VERSION 3.0)
project("SimpleMATIO" VERSION 1.1 LANGUAGES C CXX)

set(matioMinVersion 1.5.6)
set(matioDownloadVersion 1.5.11)
set(matioDownloadChecksum 41cf72778ec0a6338768f8947b88dc73)
set(matioDownloadChecksumAlgorithm MD5)
set(matioUrl "http://downloads.sourceforge.net/project/matio/matio/${matioDownloadVersion}/matio-${matioDownloadVersion}.7z")

set(zlibDownloadVersion 1.2.11)
set(zlibDownloadChecksum 4ff941449631ace0d4d203e3483be9dbc9da454084111f97ea0a2114e19bf066)
set(zlibDownloadChecksumAlgorithm SHA256)
set(zlibUrl "http://zlib.net/zlib-${zlibDownloadVersion}.tar.xz")

set(debuggingThis OFF)
option(alwaysDownloadMatio "Whether to always download libmatio regardless of any installed version" ON)
option(alwaysDownloadZlib "Whether to always download zlib regardless of any installed version" ON)




if(debuggingThis)
    function(debug_message text)
        message(${text})
    endfunction()
else()
    function(debug_message text)
    endfunction()
endif()



set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "RelWithDebugInfo")
endif()
string(TOLOWER "${CMAKE_BUILD_TYPE}" buildType)


set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)


set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)


# Do not export symbols by default
cmake_policy(SET CMP0063 NEW)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# To use bindir and others
include(GNUInstallDirs)

# where to install stuff
set(headerInstallDir "${CMAKE_INSTALL_FULL_INCLUDEDIR}/simplematio")
set(libInstallDir "${CMAKE_INSTALL_FULL_LIBDIR}/simplematio")

include(ExternalProject)

# Tell the code that we're building a library
add_definitions(-DBUILDING_SMIO)

# Where to install the pkg-config files
if(UNIX)
    set(pkgConfigInstallPath "${CMAKE_INSTALL_PREFIX}/share/pkgconfig")
else()
    set(pkgConfigInstallPath "${headerInstallDir}")
endif()

# Check dependecies
find_package(Eigen3 REQUIRED)

set(noZlib ${alwaysDownloadZlib})
if(NOT alwaysDownloadZlib)
    find_package(ZLIB)
    if(NOT ZLIB_FOUND)
        set(noZlib ON)
        message("Didn't find any Zlib, building it myself")
        message("")
    else()
        set(staticZlib "")
        foreach(lib ${ZLIB_LIBRARIES})
            get_filename_component(libExt "${lib}" EXT)
            if("${libext}" STREQUAL "${CMAKE_STATIC_LIBRARY_SUFFIX}")
                set(staticZlib "${lib}")
            endif()
        endforeach()

        if(NOT staticZlib)
            message(FATAL_ERROR "Didn't find static zlib!")
        endif()
        set(zlibName "${staticZlib}")
        list(GET ZLIB_INCLUDE_DIRS 0 zlibIncludeDir)
    endif()
endif()

if(noZlib)
    set(zlibPrefix "${CMAKE_CURRENT_BINARY_DIR}/zlib")
    set(zlibLibDir "${zlibPrefix}/lib")
    set(zlibIncludeDir "${zlibPrefix}/include")
    set(configurePrefix "${zlibPrefix}/unneeded")

    ExternalProject_Add(zlib
        URL ${zlibUrl}
        URL_HASH ${zlibDownloadChecksumAlgorithm}=${zlibDownloadChecksum}
        PREFIX ${zlibPrefix}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure
            --static --prefix="${configurePrefix}" --libdir="${zlibLibDir}" --includedir="${zlibIncludeDir}"
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make
        INSTALL_COMMAND make install
    )
    set(staticZlib "${zlibLibDir}/libz${CMAKE_STATIC_LIBRARY_SUFFIX}")


    install(
        FILES "${zlibIncludeDir}/zconf.h" "${zlibIncludeDir}/zlib.h"
        DESTINATION "${headerInstallDir}"
    )
    install(
        FILES "${staticZlib}" DESTINATION "${libInstallDir}"
    )
    set(zlibName "${libInstallDir}/libz${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif()




set(noMatio ${alwaysDownloadMatio})

if(NOT alwaysDownloadMatio)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(MATIO matio)
    if(NOT MATIO_FOUND)
        set(noMatio ON)
        message("Didn't find any libmatio, building it myself.")
        message("")
    elseif(MATIO_VERSION VERSION_LESS matioMinVersion)
        set(noMatio ON)
        message("Version ${MATIO_VERSION} of libmatio is too old, at least ${matioMinVersion} is needed. Building it myself.")
    else()
        find_library(staticLibmatio
            NAMES libmatio${CMAKE_STATIC_LIBRARY_SUFFIX}
            NO_DEFAULT_PATH
            PATHS ${MATIO_LIBRARY_DIRS}
        )
        debug_message("MATIO_LIBRARY_DIRS = ${MATIO_LIBRARY_DIRS}")
        debug_message("staticLibmatio = ${staticLibmatio}")

        if(NOT staticLibmatio)
            message(FATAL_ERROR "Didn't find static libmatio!")
        endif()
    set(libmatioName "${staticLibmatio}")
    endif()
endif()

if (noMatio)
    set(matioPrefix "${CMAKE_CURRENT_BINARY_DIR}/libmatio")
    set(matioLibDir "${matioPrefix}/lib")
    # Copy relevant zlib files to make sure that matios configure picks up the right ones
    set(matioZlibDir "${matioPrefix}/zlib")
    set(matioStaticZlib "${matioZlibDir}/lib/libz${CMAKE_STATIC_LIBRARY_SUFFIX}")
    set(matioIncludeDir "${matioPrefix}/include")
    set(configurePrefix "${matioPrefix}/unneeded")

    add_custom_target(
        copyZlibToMatioDir
        COMMAND ${CMAKE_COMMAND} -E copy "${staticZlib}" "${matioStaticZlib}"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${zlibIncludeDir}" "${matioZlibDir}/include"
        DEPENDS "${staticZlib}"
    )
    if(noZlib)
        add_dependencies(copyZlibToMatioDir zlib)
    endif()

    ExternalProject_Add(matio
        URL ${matioUrl}
        URL_HASH ${matioDownloadChecksumAlgorithm}=${matioDownloadChecksum}
        PREFIX ${matioPrefix}
        CONFIGURE_COMMAND chmod +x <SOURCE_DIR>/configure
        COMMAND <SOURCE_DIR>/configure --with-zlib=${matioZlibDir}
            --prefix=${configurePrefix} --libdir=${matioLibDir} --includedir=${matioIncludeDir} --oldincludedir=${matioIncludeDir}
            --enable-static --disable-shared --with-pic --with-default-file-ver=5
            --disable-mat73 --enable-extended-sparse --disable-profile
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make
        INSTALL_COMMAND make install
    )

    add_dependencies(matio copyZlibToMatioDir)

    install(
        FILES "${matioIncludeDir}/matio.h" "${matioIncludeDir}/matio_pubconf.h"
        DESTINATION "${headerInstallDir}"
    )
    set(staticLibmatio "${matioLibDir}/libmatio${CMAKE_STATIC_LIBRARY_SUFFIX}")
    install(FILES "${staticLibmatio}" DESTINATION "${libInstallDir}")
    set(libmatioName "${libInstallDir}/libmatio${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif()



# EIGEN3_INCLUDE_DIR without trailing S!
set(libs ${staticLibmatio} ${ZLIB_LIBRARIES})
debug_message("libs=${libs}")


# Compilation and linking directives

foreach(target simplematio simplematiof)
    add_library(${target} STATIC simplematio.c)
    set_target_properties(${target} PROPERTIES PREFIX "lib") # Always prefix with lib, even on Windows
    target_include_directories(${target} PUBLIC "${matioIncludeDir}" "${zlibIncludeDir}" "${MATIO_INCLUDE_DIRS}" "${EIGEN3_INCLUDE_DIR}" "${ZLIB_INCLUDE_DIRS}" "${CMAKE_CURRENT_SOURCE_DIR}")
    if(noMatio)
        add_dependencies(${target} matio)
    endif()
endforeach()

target_compile_definitions(simplematiof PRIVATE SMIO_USE_FLOATS=1)

foreach(target IN ITEMS simplematio simplematiof)
    target_link_libraries("${target}" "${staticLibmatio}" "${staticZlib}")
endforeach()

add_executable(EigenSimplematioTest EXCLUDE_FROM_ALL EigenSimplematioTest.cpp)
add_executable(EigenSimplematioTestFloat EXCLUDE_FROM_ALL EigenSimplematioTest.cpp)
target_compile_definitions(EigenSimplematioTestFloat PRIVATE SMIO_USE_FLOATS=1)
target_link_libraries(EigenSimplematioTest simplematio)
target_link_libraries(EigenSimplematioTestFloat simplematiof)
add_custom_target(check COMMAND true)
add_test(EigenSimplematioTest EigenSimplematioTest)
add_test(EigenSimplematioTestFloat EigenSimplematioTestFloat)
add_dependencies(check EigenSimplematioTest EigenSimplematioTestFloat)

foreach(thing "simplematio;double;;" "simplematiof;float; -DSMIO_USE_FLOATS=1;")
    list(GET thing 0 target)
    list(GET thing 1 type)
    list(GET thing 2 maybeFloatMacro)
    set(configuredFile "${CMAKE_CURRENT_BINARY_DIR}/${target}.pc")

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pkgConfigTemplate.pc" "${configuredFile}" @ONLY NEWLINE_STYLE UNIX)
    install(FILES "${configuredFile}"
        DESTINATION "${pkgConfigInstallPath}"
    )
endforeach()


install(TARGETS simplematio simplematiof
    ARCHIVE DESTINATION "${libInstallDir}"
)
install(FILES simplematio.h EigenSimplematio.hpp EigenSimplematio.cpp
    DESTINATION "${headerInstallDir}")
