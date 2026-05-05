include (GoogleTest)

set (TestNum 0)

# test executables
macro (ADD_TEST test_file)
    # set test name
    string (REGEX MATCH "([^\/]+)\.cpp" _ "${test_file}")
    set (test_name "${TestNum}_${CMAKE_MATCH_1}")
    math (EXPR TestNum "${TestNum}+1")

    maki_config_print ("Adding test '${test_name}' @ ${CMAKE_CURRENT_LIST_DIR}/${test_file}.")
    
    # create test target
    add_executable (${test_name} ${MAKI_CLASSIFY_SOURCES} "${CMAKE_CURRENT_LIST_DIR}/${test_file}")
    target_link_libraries (${test_name} maki::maki)
    gtest_discover_tests (${test_name})
    list (APPEND MAKI_EXECUTABLES ${test_name})
endmacro ()

# test suite

ADD_TEST ("build/seqio/test_detect_file_types.cpp")
ADD_TEST ("build/seqio/test_gff_parser.cpp")
ADD_TEST ("build/seqio/test_gff_iter.cpp")
ADD_TEST ("build/seqio/test_chunked_fasta.cpp")

ADD_TEST ("build/colours/test_buffer_value.cpp")
ADD_TEST ("build/colours/test_colours.cpp")
ADD_TEST ("build/colours/test_metacolours.cpp")
ADD_TEST ("build/colours/test_construct_registry.cpp")

ADD_TEST ("build/suffixes/test_suffix_encoding.cpp")
ADD_TEST ("build/suffixes/test_suffix_counting.cpp")
ADD_TEST ("build/suffixes/parallel/test_suffix_gate.cpp")
ADD_TEST ("build/suffixes/test_kmer_buffer.cpp")
ADD_TEST ("build/suffixes/test_terminal_buffer.cpp")
ADD_TEST ("build/suffixes/test_sorting.cpp")

ADD_TEST ("build/archives/test_packed_array.cpp")
ADD_TEST ("build/archives/test_archive_basic.cpp")
ADD_TEST ("build/archives/test_archive_validation.cpp")
ADD_TEST ("build/archives/test_sdsl_sinks.cpp")
ADD_TEST ("build/archives/test_bytearray_sink.cpp")
ADD_TEST ("build/archives/test_multisink_pipeline.cpp")

ADD_TEST ("build/wm/test_wm_construct.cpp")
ADD_TEST ("build/graph/test_push.cpp")
ADD_TEST ("build/graph/test_construct.cpp")

ADD_TEST ("core/graph/test_navigate.cpp")

ADD_TEST ("classify/seqio/test_fastq_parser.cpp")
ADD_TEST ("classify/interleave/test_elm_merge_small.cpp")
ADD_TEST ("classify/interleave/test_elm_merge_large.cpp")
ADD_TEST ("classify/interleave/test_merge.cpp")

ADD_TEST ("core/utils/test_wap_vector.cpp")
ADD_TEST ("core/graph/archive/test_counts.cpp")

## trees

#ADD_TEST ("tree/test_newick.cpp")

# append executables in parent scope
set (MAKI_EXECUTABLES ${MAKI_EXECUTABLES} PARENT_SCOPE)
