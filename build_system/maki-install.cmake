# Sean Solari - sean.solari@monash.edu
# Microbiota and Systems Biology Lab, CIIID, Hudson Institute of Medical Research

cmake_minimum_required(VERSION 3.21)

include (GNUInstallDirs)

# install cmake files in /share/cmake
install (FILES "${MAKI_CLONE_DIR}/build_system/maki-config.cmake"
               "${MAKI_CLONE_DIR}/build_system/maki-config-version.cmake"
         DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/maki")

# install seqan3 header files in /include/seqan3
install (DIRECTORY "${MAKI_INCLUDE_DIR}/maki" TYPE INCLUDE)
