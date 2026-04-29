
set (MAKI_BUILD_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/graph/wm/wm_construct.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/graph/build_colours.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/graph/construct_common.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/graph/construct_cdbg.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/graph/construct_wdbg.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/graph/interleave_buffers.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/io/fasta.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/kmers/buffers/kmers.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/kmers/buffers/nt_encoding.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/kmers/buffers/sort.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/kmers/buffers/terminals.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/kmers/construct_terminals.cpp"
  PARENT_SCOPE
)
