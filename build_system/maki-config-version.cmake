# Sean Solari - sean.solari@monash.edu
# Microbiota and Systems Biology Lab, CIIID, Hudson Institute of Medical Research

# find clone and include directories
find_path (MAKI_CLONE_DIR
           NAMES build_system/maki-config.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/..")
find_path (MAKI_INCLUDE_DIR
           NAMES maki/version.hpp
           HINTS "${MAKI_CLONE_DIR}/include")
# external header-only libs
find_path (MAKI_DEPENDENCY_INCLUDE_DIRS
           NAMES sdsl/version.hpp
           HINTS "${MAKI_CLONE_DIR}/extern")

# extract version from maki/version.hpp header
file (STRINGS "${MAKI_INCLUDE_DIR}/maki/version.hpp" MAKI_VERSION_HPP
      REGEX "#define MAKI_VERSION_(MAJOR|MINOR|PATCH)")
string (REGEX REPLACE "#define MAKI_VERSION_(MAJOR|MINOR|PATCH) " "" PACKAGE_VERSION "${MAKI_VERSION_HPP}")
string (REGEX REPLACE ";" "." PACKAGE_VERSION "${PACKAGE_VERSION}")

if (PACKAGE_VERSION VERSION_LESS PACKAGE_FIND_VERSION)
    set (PACKAGE_VERSION_COMPATIBLE FALSE)
else ()

    if (PACKAGE_VERSION MATCHES "^([0-9]+)\\.")
        set (_PACKAGE_VERSION_MAJOR "${CMAKE_MATCH_1}")
    endif ()

    if (PACKAGE_FIND_VERSION_MAJOR STREQUAL _PACKAGE_VERSION_MAJOR)
        set (PACKAGE_VERSION_COMPATIBLE TRUE)
    else ()
        set (PACKAGE_VERSION_COMPATIBLE FALSE)
    endif ()

    if (PACKAGE_FIND_VERSION STREQUAL PACKAGE_VERSION)
        set (PACKAGE_VERSION_EXACT TRUE)
    endif ()
endif ()

# extract release candidate
file (STRINGS "${MAKI_INCLUDE_DIR}/maki/version.hpp" MAKI_RELEASE_CANDIDATE_HPP
      REGEX "#define MAKI_RELEASE_CANDIDATE ")
string (REGEX REPLACE "#define MAKI_RELEASE_CANDIDATE " "" MAKI_RELEASE_CANDIDATE_VERSION
                      "${MAKI_RELEASE_CANDIDATE_HPP}")

# MAKI_PROJECT_VERSION is intended to be used within `project (... VERSION "${MAKI_PROJECT_VERSION}")`.
set (MAKI_PROJECT_VERSION "${PACKAGE_VERSION}")
if (MAKI_RELEASE_CANDIDATE_VERSION VERSION_GREATER "0")
    set (PACKAGE_VERSION "${PACKAGE_VERSION}-rc.${MAKI_RELEASE_CANDIDATE_VERSION}")
endif ()

if (NOT MAKI_PROJECT_VERSION VERSION_EQUAL PACKAGE_VERSION)
    message (AUTHOR_WARNING "MAKI_PROJECT_VERSION and MAKI_VERSION mismatch, "
                            "please report this issue and mention your cmake version.")
endif ()
