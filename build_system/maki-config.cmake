# Sean Solari - sean.solari@monash.edu
# Microbiota and Systems Biology Lab, CIIID, Hudson Institute of Medical Research
#
# This CMake module will try to find maki and its dependencies.
#
# Usage:
#   find_package (maki [REQUIRED] ...)
#
# maki has the following platform requirements:
#
#   C++20
#
# maki has the following library requirements:
#
#   SDSL            - Succinct Data Structures Library
#   cereal          - Serialisation of succinct data structures
#   GLog            - Logging
#   GTest           - Unit Test construction suite
#   oneapi::tbb     - Intel Threading Building Blocks
#   SeqAn3          - Sequence Analysis C++ Toolkit
#   
# Once maki has been found, the following variables will be set.
#
#   MAKI_FOUND
#
#   MAKI_VERSION
#   MAKI_VERSION_MAJOR
#   MAKI_VERSION_MINOR
#   MAKI_VERSION_PATCH
#
#   MAKI_INCLUDE_DIRS
#   MAKI_LIBRARIES
#   MAKI_DEFINITIONS
#   MAKI_CXX_FLAGS
#
# The following targets are defined.
#
#   maki::maki      -- interface target where
#                               target_link_libraries([target] maki::maki)
#                          automatically sets
#                               target_include_directories([target] $MAKI_INCLUDE_DIRS),
#                               target_link_libraries([target] $MAKI_LIBRARIES),
#                               target_compile_definitions([target] $MAKI_DEFINITIONS) and
#                               target_compile_options([target] $MAKI_CXX_FLAGS)
#
# ================================================================================

cmake_minimum_required(VERSION 3.21)

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

# make output globally quiet if required by find_package, this effects cmake functions like `check_*`
set (CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set (CMAKE_REQUIRED_QUIET ${${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY})

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Finding maki and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (ExternalProject)
include (FetchContent)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (maki_config_print text)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (maki_config_error text)
    if (${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
        message (FATAL_ERROR ${text})
    else ()
        if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message (WARNING ${text})
        endif ()
        return ()
    endif ()
endmacro ()

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}Finding maki and checking requirements:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Find maki include path
# ----------------------------------------------------------------------------

# get path to compiler
maki_config_print ("Using compiler ${CMAKE_CXX_COMPILER}")
get_filename_component (_CMAKE_CXX_COMPILER_BIN ${CMAKE_CXX_COMPILER} DIRECTORY)

# append to prefix search path
list (APPEND CMAKE_PREFIX_PATH ${_CMAKE_CXX_COMPILER_BIN}/../)
maki_config_print ("CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")

# * `MAKI_INCLUDE_DIR` was already found in maki-config-version.cmake
if (MAKI_INCLUDE_DIR)
    maki_config_print ("maki include dir found:   ${MAKI_INCLUDE_DIR}")
else ()
    maki_config_error ("maki include directory could not be found (MAKI_INCLUDE_DIR: '${MAKI_INCLUDE_DIR}')")
endif ()
if (MAKI_DEPENDENCY_INCLUDE_DIRS)
    maki_config_print ("maki dependencies include dir found:   ${MAKI_DEPENDENCY_INCLUDE_DIRS}")
else ()
    maki_config_error ("maki dependencies include directory could not be found (MAKI_DEPENDENCY_INCLUDE_DIRS: '${MAKI_DEPENDENCY_INCLUDE_DIRS}')")
endif ()

# ----------------------------------------------------------------------------
# Collect submodules
# ----------------------------------------------------------------------------

# use global variables in Check* calls
list (APPEND CMAKE_INCLUDE_PATH ${MAKI_INCLUDE_DIR} ${MAKI_DEPENDENCY_INCLUDE_DIRS})
set (CMAKE_REQUIRED_INCLUDES ${CMAKE_INCLUDE_PATH} )
set (CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS})

# ----------------------------------------------------------------------------
# Require C++20
# ----------------------------------------------------------------------------

set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# ----------------------------------------------------------------------------
# Compiler Compatability
# ----------------------------------------------------------------------------

set (gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set (msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")

# ----------------------------------------------------------------------------
# Declare all dependencies for future FetchContent_MakeAvailable
# ----------------------------------------------------------------------------

# SeqAn3 --

FetchContent_Declare (
    seqan3_fetch_content
    GIT_REPOSITORY https://github.com/seqan/seqan3.git
    GIT_TAG        c8e7fdbe2aa259b96d0007ece0218d0150fca045 # v3.4.0-rc.3
    )

# gtl --

FetchContent_Declare(
    gtl
    GIT_REPOSITORY https://github.com/greg7mdp/gtl.git
    GIT_TAG        v1.2.0 # adjust tag/branch/commit as needed
    )

# ZStr --

FetchContent_Declare(
    ZStrGitRepo
    GIT_REPOSITORY    "https://github.com/mateidavid/zstr"
    GIT_TAG           "master"
    )

# GTest --

FetchContent_Declare (
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        f8d7d77c06936315286eb55f8de22cd23c188571 # release-1.14.0
    FIND_PACKAGE_ARGS NAMES GTest
    )

# For Windows: Prevent overriding the parent project's compiler/linker settings
set (gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# GLog --

FetchContent_Declare (
    googlelogging
    GIT_REPOSITORY https://github.com/google/glog.git
    GIT_TAG        34b8da6496aec6a98277808701cfa834fae9801f # release-0.7.0
    FIND_PACKAGE_ARGS NAMES glog
    )

# TBB --

FetchContent_Declare (
    onetbb
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG        45587e94dfb6dfe00220c5f520020a5bc745e92f # release-2022.1.0
    FIND_PACKAGE_ARGS NAMES TBB
    )

set (TBB_STRICT OFF)

##

FetchContent_MakeAvailable (
    seqan3_fetch_content
    gtl
    ZStrGitRepo
    googletest
    onetbb
    )

# ----------------------------------------------------------------------------
# Require SDSL - Succinct Data Structures Library
# ----------------------------------------------------------------------------

check_include_file_cxx (sdsl/version.hpp _MAKI_HAVE_SDSL)

if (_MAKI_HAVE_SDSL)
    maki_config_print ("Required dependency:        SDSL found.")
else ()
    maki_config_error (
        "The SDSL library is required, but wasn't found. Get it from https://github.com/xxsds/sdsl-lite")
endif ()

# ----------------------------------------------------------------------------
# Require Cereal - Serialisation of succinct data structures
# ----------------------------------------------------------------------------

# Note the cereal::cereal library defined in cereal/CMakeLists.txt is encompassed
# by including all needed headers, which is checked here.
check_include_file_cxx (cereal/version.hpp _MAKI_HAVE_CEREAL)

if (_MAKI_HAVE_CEREAL)
    maki_config_print ("Required dependency:        Cereal found.")
    set (MAKI_DEFINITIONS SDSL_HAS_CEREAL ${MAKI_DEFINITIONS})
else ()
    maki_config_error (
        "The Cereal library is required, but wasn't found. Get it from https://github.com/USCiLab/cereal")
endif ()

# ----------------------------------------------------------------------------
# Require GTL - Greg's Template Library of useful classes
# ----------------------------------------------------------------------------

if (TARGET gtl)
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} gtl)
    # collect desired attributes
    get_target_property (gtl_VERSION gtl VERSION)
    maki_config_print ("Required dependency:        gtl-${gtl_VERSION} found.")
else ()
    maki_config_error ("Dependency gtl not found.")
endif ()

# ----------------------------------------------------------------------------
# Require Indicators - Thread-safe progress bars and spinners
# ----------------------------------------------------------------------------

check_include_file_cxx (indicators/progress_bar.hpp _MAKI_HAVE_INDICATORS)

if (_MAKI_HAVE_INDICATORS)
    maki_config_print ("Required dependency:        Indicators found.")
else ()
    maki_config_error (
        "The Indicators library is required, but wasn't found. Get it from https://github.com/p-ranav/indicators")
endif ()

# ----------------------------------------------------------------------------
# Require SeqAn3 - Sequence Analysis C++ Toolkit
# ----------------------------------------------------------------------------

if (TARGET seqan3::seqan3)
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} seqan3::seqan3)
    maki_config_print ("Required dependency:        SeqAn3-${SEQAN3_VERSION} found.")
    maki_config_print ("SeqAn3 include dirs: ${SEQAN3_INCLUDE_DIRS}")
    maki_config_print ("SeqAn3 libs: ${SEQAN3_LIBRARIES}")
else ()
    maki_config_error ("Dependency SeqAn3 not found.")
endif ()

# ----------------------------------------------------------------------------
# Require OMP
# ----------------------------------------------------------------------------

find_package(OpenMP REQUIRED)

if (TARGET OpenMP::OpenMP_CXX)
    maki_config_print ("Required dependency:        OMP found.")
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} OpenMP::OpenMP_CXX)
else ()
    maki_config_error ("Dependency OpenMP not found.")
endif ()

# ----------------------------------------------------------------------------
# Require ZStr - C++ ZLib wrapper
# ----------------------------------------------------------------------------

find_package(ZLIB 1.2.3 REQUIRED)

if (TARGET ZLIB::ZLIB)
    maki_config_print ("Required dependency:        zlib-${ZLIB_VERSION_STRING} found.")
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} ZLIB::ZLIB)
else ()
    maki_config_error ("Dependency zlib not found.")
endif ()

if (TARGET zstr)
    # collect desired attributes
    get_target_property (zstr_VERSION zstr VERSION)
    maki_config_print ("Required dependency:        zstr-${zstr_VERSION} found.")
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} zstr::zstr)
else ()
    maki_config_error ("Dependency zstr not found.")
endif ()

# ----------------------------------------------------------------------------
# Require GTest - Unit Test construction suite
# ----------------------------------------------------------------------------

if (TARGET GTest::gtest_main)
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} GTest::gtest_main GTest::gmock_main)
    # collect desired attributes
    get_target_property (GTest_VERSION GTest::gtest_main VERSION)
    maki_config_print ("Required dependency:        googletest-${GTest_VERSION} found.")
else ()
    maki_config_error ("Dependency GTest not found.")
endif ()

# ----------------------------------------------------------------------------
# Require libatomic
# ----------------------------------------------------------------------------

find_library(ATOMIC_LIBRARY NAMES atomic atomic.so.1 libatomic libatomic.so.1)

# Check if libatomic was found
if (ATOMIC_LIBRARY)
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} ${ATOMIC_LIBRARY})
    maki_config_print ("Required dependency:        libatomic found.")
else ()
    maki_config_error ("Dependency libatomic not found.")
endif()

# ----------------------------------------------------------------------------
# Require TBB - Intel Threading Building Blocks
# ----------------------------------------------------------------------------

if (TARGET TBB::tbb)
    set (MAKI_LIBRARIES ${MAKI_LIBRARIES} TBB::tbb)
    # collect desired attributes
    get_target_property (TBB_VERSION TBB::tbb VERSION)
    maki_config_print ("Required dependency:        oneTBB-${TBB_VERSION} found.")
else ()
    maki_config_error ("Dependency TBB not found.")
endif ()

# bug fix for gcc-12 (see https://github.com/oneapi-src/oneTBB/issues/1508)
if (UNIX)
    set(TBB_WARNING_SUPPRESS ${TBB_WARNING_SUPPRESS} -Wno-stringop-overflow)
endif()

# ----------------------------------------------------------------------------
# Compile time options and build parameters
# ----------------------------------------------------------------------------

option (DO_TIMING "Activate timing methods" OFF)
option (CHECK_DATA_RACES "Build with -fsanitize=thread" OFF)
option (CHECK_OOB "Build with -fsanitize=address" OFF)

# build parameters
set (TREE_ONE_PATH       "${PROJECT_SOURCE_DIR}/data/one.nwk")
set (TREE_TWO_PATH       "${PROJECT_SOURCE_DIR}/data/two.nwk")
set (SEQUENCE_FNA_PATH   "${PROJECT_SOURCE_DIR}/data/sequence.fna")
set (SMALL_FNA_PATH      "${PROJECT_SOURCE_DIR}/data/small.fna")
set (GRAPH_FNA_PATH      "${PROJECT_SOURCE_DIR}/data/1902.02889v5.fna")
set (CONTIG_ONE_PATH     "${PROJECT_SOURCE_DIR}/data/contig_one.fna")
set (CONTIG_TWO_PATH     "${PROJECT_SOURCE_DIR}/data/contig_two.fna")
set (CONTIG_THREE_PATH   "${PROJECT_SOURCE_DIR}/data/contig_three.fna")
set (BENCHMARK_FASTA     "${PROJECT_SOURCE_DIR}/data/benchmark.fna")
set (GCF_000005845_FNA   "${PROJECT_SOURCE_DIR}/data/GCF_000005845.2_ASM584v2_genomic.fna")
set (GCF_000006765_FNA   "${PROJECT_SOURCE_DIR}/data/GCF_000006765.1_ASM676v1_genomic.fna")
set (GCF_000006925_FNA   "${PROJECT_SOURCE_DIR}/data/GCF_000006925.2_ASM692v2_genomic.fna")
set (GCF_000006945_FNA   "${PROJECT_SOURCE_DIR}/data/GCF_000006945.2_ASM694v2_genomic.fna")
set (GCF_000007765_FNA   "${PROJECT_SOURCE_DIR}/data/GCF_000007765.2_ASM776v2_genomic.fna")
set (GCF_000008725_FNA   "${PROJECT_SOURCE_DIR}/data/GCF_000008725.1_ASM872v1_genomic.fna")
set (GCA_011523445_GFF   "${PROJECT_SOURCE_DIR}/data/GCA_011523445.1_genomic.gff3")
set (GCA_021324415_GFF   "${PROJECT_SOURCE_DIR}/data/GCA_021324415.1_genomic.gff3")
set (GCA_026393455_GFF   "${PROJECT_SOURCE_DIR}/data/GCA_026393455.1_genomic.gff3")
set (GCA_026395095_GFF   "${PROJECT_SOURCE_DIR}/data/GCA_026395095.1_genomic.gff3")
set (GCA_031265645_GFF   "${PROJECT_SOURCE_DIR}/data/GCA_031265645.1_genomic.gff3")
set (GCF_000620805_GFF   "${PROJECT_SOURCE_DIR}/data/GCF_000620805.1_genomic.gff3")
set (GCF_002214165_GFF   "${PROJECT_SOURCE_DIR}/data/GCF_002214165.1_genomic.gff3")
set (GCF_003286415_GFF   "${PROJECT_SOURCE_DIR}/data/GCF_003286415.2_genomic.gff3")
set (GCF_013408985_GFF   "${PROJECT_SOURCE_DIR}/data/GCF_013408985.1_genomic.gff3")
set (GCF_020054945_GFF   "${PROJECT_SOURCE_DIR}/data/GCF_020054945.1_genomic.gff3")
set (SEQA_FNA            "${PROJECT_SOURCE_DIR}/data/seqA.fna")
set (SEQB_FNA            "${PROJECT_SOURCE_DIR}/data/seqB.fna")
set (SEQC_FNA            "${PROJECT_SOURCE_DIR}/data/seqC.fna")
set (SEQA_GFF            "${PROJECT_SOURCE_DIR}/data/seqA.gff3")
set (SEQB_GFF            "${PROJECT_SOURCE_DIR}/data/seqB.gff3")
set (SEQC_GFF            "${PROJECT_SOURCE_DIR}/data/seqC.gff3")
set (TEST_ANNOTS         "${PROJECT_SOURCE_DIR}/data/testAnnots.gff3")
set (GCF_000005845_FWD   "${PROJECT_SOURCE_DIR}/data/fq/GCF_000005845.2.1.fq.gz")
set (GCF_000005845_REV   "${PROJECT_SOURCE_DIR}/data/fq/GCF_000005845.2.2.fq.gz")
set (GCF_000006765_FWD   "${PROJECT_SOURCE_DIR}/data/fq/GCF_000006765.1.1.fq.gz")
set (GCF_000006765_REV   "${PROJECT_SOURCE_DIR}/data/fq/GCF_000006765.1.2.fq.gz")
set (GCF_000006925_FWD   "${PROJECT_SOURCE_DIR}/data/fq/GCF_000006925.2.1.fq.gz")
set (GCF_000006925_REV   "${PROJECT_SOURCE_DIR}/data/fq/GCF_000006925.2.2.fq.gz")
set (GCF_000005845_SMFWD "${PROJECT_SOURCE_DIR}/data/fq/small/GCF_000005845.2.1.fq")
set (GCF_000005845_SMREV "${PROJECT_SOURCE_DIR}/data/fq/small/GCF_000005845.2.2.fq")

# data type parameters
set (BITS_PER_EDGE               3)
set (KMER_OVERLAP_GRAINSIZE      100000)
set (KMER_RADIX_CHUNKSIZE        100000)
set (PARALLEL_FLUSH_TOKENS       100)
set (GRAPH_FLUSH_MAX_TEMP_BYTES  100000000)
set (ZSTD_COMPRESSION_LEVEL      3)

# ----------------------------------------------------------------------------
# Compiler flags
# ----------------------------------------------------------------------------

# specify the C++ standard and compiler flags
add_library (maki_compiler_flags INTERFACE)
target_compile_features (maki_compiler_flags INTERFACE cxx_std_20)

# project-wide settings
set_target_properties (maki_compiler_flags PROPERTIES
                       INTERFACE_EXPORT_COMPILE_COMMANDS ON
                       INTERFACE_VERBOSE_MAKEFILE ON
                       INTERFACE_CXX_STANDARD_REQUIRED ON
                       INTERFACE_CXX_EXTENSIONS OFF
                       )

# add compiler warning flags when building
target_compile_options (maki_compiler_flags INTERFACE
                        "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wformat=2;-Wunused>>"
                        "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
                        "-Wno-interference-size"
                        )

# check threading
if (CHECK_DATA_RACES)
    target_compile_options (maki_compiler_flags INTERFACE "-fsanitize=thread")
    target_link_options (maki_compiler_flags INTERFACE "-fsanitize=thread")
endif ()

# check out of bounds access/write
if (CHECK_OOB)
    target_compile_options (maki_compiler_flags INTERFACE "-fsanitize=address")
    target_link_options (maki_compiler_flags INTERFACE "-fsanitize=address")
endif ()

# check timing
if (DO_TIMING)
    set (USE_TIMING_ROUTINES 1)
else ()
    set (USE_TIMING_ROUTINES 0)
endif ()

maki_config_print ("Using timing routines?    ${USE_TIMING_ROUTINES} (bool)")

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS MAKI_INCLUDE_DIR)

# Set MAKI_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
# This needs to be done, because `find_package(maki)` might be called in any case-sensitive way and we want to
# guarantee that MAKI_* are always set.
foreach (package_var
         FOUND
         DIR
         ROOT
         CONFIG
         VERSION
         VERSION_MAJOR
         VERSION_MINOR
         VERSION_PATCH
         VERSION_TWEAK
         VERSION_COUNT)
    set (MAKI_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
endforeach ()

# propagate MAKI_INCLUDE_DIR into MAKI_INCLUDE_DIRS
set (MAKI_INCLUDE_DIRS ${MAKI_INCLUDE_DIR} ${MAKI_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

if (MAKI_FOUND AND NOT TARGET maki::maki)
    maki_config_print ("maki include dir: ${MAKI_INCLUDE_DIR} ")
    maki_config_print ("maki dependency include dirs: ${MAKI_DEPENDENCY_INCLUDE_DIRS}")
    separate_arguments (MAKI_CXX_FLAGS_LIST UNIX_COMMAND "${MAKI_CXX_FLAGS}")

    add_library (maki_maki INTERFACE)
    target_include_directories (maki_maki INTERFACE
                                "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>" # needed for configured 'maki.h'
                                "$<INSTALL_INTERFACE:include>"
                                ${MAKI_INCLUDE_DIR}
                                ${MAKI_DEPENDENCY_INCLUDE_DIRS}
                                )
    target_link_libraries (maki_maki INTERFACE maki_compiler_flags ${MAKI_LIBRARIES})
    target_compile_definitions (maki_maki INTERFACE ${MAKI_DEFINITIONS})
    target_compile_options (maki_maki INTERFACE ${MAKI_CXX_FLAGS_LIST})
    add_library (maki::maki ALIAS maki_maki)
endif ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (MAKI_FIND_DEBUG)
    message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
    message ("")
    message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
    message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
    message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
    message ("  MAKI_INCLUDE_DIR          ${MAKI_INCLUDE_DIR}")
    message ("")
    message ("")
    message ("  MAKI_INCLUDE_DIRS         ${MAKI_INCLUDE_DIRS}")
    message ("  MAKI_LIBRARIES            ${MAKI_LIBRARIES}")
    message ("  MAKI_DEFINITIONS          ${MAKI_DEFINITIONS}")
    message ("  MAKI_CXX_FLAGS            ${MAKI_CXX_FLAGS}")
    message ("")
    message ("  MAKI_VERSION              ${MAKI_VERSION}")
    message ("  MAKI_VERSION_MAJOR        ${MAKI_VERSION_MAJOR}")
    message ("  MAKI_VERSION_MINOR        ${MAKI_VERSION_MINOR}")
    message ("  MAKI_VERSION_PATCH        ${MAKI_VERSION_PATCH}")
endif ()
