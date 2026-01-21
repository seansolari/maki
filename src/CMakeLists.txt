
add_subdirectory ("core")
add_subdirectory ("build")
add_subdirectory ("classify")

list (APPEND MAKI_BUILD_SOURCES ${MAKI_CORE_SOURCES})
list (APPEND MAKI_CLASSIFY_SOURCES ${MAKI_BUILD_SOURCES})
