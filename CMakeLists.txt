# This file provides capactiy to download the dependencies of maki, build the
# main execution binaries, as well as its test modules.

cmake_minimum_required (VERSION 3.24)

find_path (MAKI_MODULE_PATH "maki-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/build_system")
list (APPEND CMAKE_MODULE_PATH "${MAKI_MODULE_PATH}" "${PROJECT_SOURCE_DIR}/cmake")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_VERBOSE_MAKEFILE 0)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

include (maki-config-version)

project (maki
         LANGUAGES CXX
         VERSION 2.0.0
         DESCRIPTION "metagenomics using an annotated k-mer index"
         HOMEPAGE_URL "https://github.com/seansolari/maki-engine"
         )

find_package (maki 2.0 REQUIRED HINTS ${MAKI_MODULE_PATH})

maki_config_print ("Using maki-${MAKI_VERSION}.")

# configurable options
option (ENABLE_MAKI_TESTING "Enable testing suite" OFF)

# configurable header
configure_file (cmake/maki.h.in include/maki/maki.h)
configure_file (${CMAKE_SOURCE_DIR}/cmake/CTestCustom.cmake ${CMAKE_BINARY_DIR} @ONLY)

# collect sources
add_subdirectory ("src")

# custom test targets
if (ENABLE_MAKI_TESTING)
  add_subdirectory ("test")

  # custom build target to build maki-related targets
  add_custom_target (maki_build_targets)
  add_dependencies (maki_build_targets ${MAKI_EXECUTABLES})
endif ()

# Python bindings

set (PYBIND11_FINDPYTHON ON)
find_package (pybind11 REQUIRED)

pybind11_add_module (_maki ${MAKI_CLASSIFY_SOURCES} "${MAKI_CLONE_DIR}/bindings/python/maki.cpp")
target_link_libraries (_maki PRIVATE maki::maki)

install (TARGETS _maki
  LIBRARY DESTINATION maki
)
