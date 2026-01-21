# Sean Solari - sean.solari@monash.edu
# Microbiota and Systems Biology Lab, CIIID, Hudson Institute of Medical Research

# This file provides capactiy to download the dependencies of maki, build the
# main execution binaries, as well as its test modules.

cmake_minimum_required (VERSION 3.24)

find_path (MAKI_MODULE_PATH "maki-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/build_system")
list (APPEND CMAKE_MODULE_PATH "${MAKI_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_VERBOSE_MAKEFILE 0)

include (maki-config-version)

project (maki
         LANGUAGES CXX
         VERSION 2.0.0
         DESCRIPTION "microbial annotated k-mer index"
         HOMEPAGE_URL "https://github.com/seansolari/pam2s"
         )

find_package (maki 2.0 REQUIRED HINTS ${MAKI_MODULE_PATH})

maki_config_print ("Using maki-${MAKI_VERSION}.")

# configurable options
option (INSTALL_MAKI "Enable installation of maki binaries." OFF)
option (ENABLE_MAKI_TESTING "Enable testing suite" ON)

# configurable header
configure_file (cmake/maki.h.in include/maki/maki.h)
configure_file (${CMAKE_SOURCE_DIR}/cmake/CTestCustom.cmake ${CMAKE_BINARY_DIR} @ONLY)

# collect sources
add_subdirectory ("src")

# two executables
# 1) maki-build
add_executable(maki-build ${MAKI_BUILD_SOURCES} "${MAKI_CLONE_DIR}/src/build.cpp")
target_link_libraries (maki-build maki::maki)

# 2) maki-clsasify
add_executable(maki-classify ${MAKI_CLASSIFY_SOURCES} "${MAKI_CLONE_DIR}/src/classify.cpp")
target_link_libraries (maki-classify maki::maki)

list (APPEND MAKI_EXECUTABLES maki-build maki-classify)

# custom test targets
if (ENABLE_MAKI_TESTING)
    add_subdirectory ("test/maki")
endif ()

# custom build target to build maki-related targets
add_custom_target (maki_build_targets)
add_dependencies (maki_build_targets ${MAKI_EXECUTABLES})

# intall binaries, headers and libs
if (INSTALL_MAKI)
    include (maki-install)
endif ()
