#pragma once
// Sean Solari - sean.solari@monash.edu, 2024
// Microbiota and Systems Biology Lab, CIIID, Hudson Institute of Medical Research

#include <cstddef>
#include <cstdint>

//!\brief The major version as MACRO.
#define MAKI_VERSION_MAJOR 2
//!\brief The minor version as MACRO.
#define MAKI_VERSION_MINOR 0
//!\brief The patch version as MACRO.
#define MAKI_VERSION_PATCH 1
//!\brief The release candidate number. 0 means stable release, >= 1 means release candidate.
#define MAKI_RELEASE_CANDIDATE 1

//!\brief The full version as MACRO (number).
#define MAKI_VERSION (MAKI_VERSION_MAJOR * 10000 + MAKI_VERSION_MINOR * 100 + MAKI_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define MAKI_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define MAKI_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH) \
    MAKI_VERSION_CSTRING_HELPER_STR(MAJOR)                    \
    "." MAKI_VERSION_CSTRING_HELPER_STR(MINOR) "." MAKI_VERSION_CSTRING_HELPER_STR(PATCH)

#if (MAKI_RELEASE_CANDIDATE > 0)
//!\brief A helper function that expands to a suitable release candidate suffix.
#define MAKI_RELEASE_CANDIDATE_HELPER(RC) "-rc." MAKI_VERSION_CSTRING_HELPER_STR(RC)
#else
//!\brief A helper function that expands to a suitable release candidate suffix.
#define MAKI_RELEASE_CANDIDATE_HELPER(RC) ""
#endif

//!\brief The full version as null terminated string.
#define MAKI_VERSION_CSTRING                                                                           \
    MAKI_VERSION_CSTRING_HELPER_FUNC(MAKI_VERSION_MAJOR, MAKI_VERSION_MINOR, MAKI_VERSION_PATCH) \
    MAKI_RELEASE_CANDIDATE_HELPER(MAKI_RELEASE_CANDIDATE)

namespace MAKI
{

    //!\brief The major version.
    constexpr uint8_t MAKI_version_major = MAKI_VERSION_MAJOR;
    //!\brief The minor version.
    constexpr uint8_t MAKI_version_minor = MAKI_VERSION_MINOR;
    //!\brief The patch version.
    constexpr uint8_t MAKI_version_patch = MAKI_VERSION_PATCH;

    //!\brief The full version as `std::size_t`.
    constexpr std::size_t MAKI_version = MAKI_VERSION;

    //!\brief The full version as null terminated string.
    constexpr char const *MAKI_version_cstring = MAKI_VERSION_CSTRING;

} // namespace maki
