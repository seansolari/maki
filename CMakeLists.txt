# Sean Solari - sean.solari@monash.edu
# Microbiota and Systems Biology Lab, CIIID, Hudson Institute of Medical Research

# This file provides capactiy to download the dependencies of maki, build the
# main execution binaries, as well as its test modules.

cmake_minimum_required (VERSION 3.24)

find_path (MAKI_MODULE_PATH "maki-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/build_system")
list (APPEND CMAKE_MODULE_PATH "${MAKI_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_VERBOSE_MAKEFILE 0)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
option (ENABLE_MAKI_TESTING "Enable testing suite" ON)

# configurable header
configure_file (cmake/maki.h.in include/maki/maki.h)
configure_file (${CMAKE_SOURCE_DIR}/cmake/CTestCustom.cmake ${CMAKE_BINARY_DIR} @ONLY)

# collect sources
add_subdirectory ("src")

message (STATUS "`Build` sources: ${MAKI_BUILD_SOURCES}")
message (STATUS "`Classify` sources: ${MAKI_CLASSIFY_SOURCES}")

# custom test targets
if (ENABLE_MAKI_TESTING)
    add_subdirectory ("test")
endif ()

# custom build target to build maki-related targets
add_custom_target (maki_build_targets)
add_dependencies (maki_build_targets ${MAKI_EXECUTABLES} )

# Python bindings

set (PYBIND11_FINDPYTHON ON)
find_package (pybind11 REQUIRED)

pybind11_add_module (_maki ${MAKI_CLASSIFY_SOURCES} "${MAKI_CLONE_DIR}/bindings/python/maki.cpp")
target_link_libraries (_maki PRIVATE maki::maki)

set_target_properties(_maki PROPERTIES
  BUILD_RPATH "$ORIGIN/libs"
  INSTALL_RPATH "$ORIGIN/libs"
)

install(TARGETS _maki
    LIBRARY DESTINATION maki
)
