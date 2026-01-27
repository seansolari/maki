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

ADD_TEST ("build/colours/test_buffer_value.cpp")
ADD_TEST ("build/colours/test_colours.cpp")
ADD_TEST ("build/colours/test_metacolours.cpp")
ADD_TEST ("build/colours/test_construct_registry.cpp")

ADD_TEST ("build/archives/test_packed_array.cpp")
ADD_TEST ("build/archives/test_archive_basic.cpp")
ADD_TEST ("build/archives/test_archive_validation.cpp")
ADD_TEST ("build/archives/test_sdsl_sinks.cpp")
ADD_TEST ("build/archives/test_bytearray_sink.cpp")
ADD_TEST ("build/archives/test_multisink_pipeline.cpp")

ADD_TEST ("build/wm/test_wm_construct.cpp")

# parsing

#ADD_TEST ("io/test_parsing.cpp")
#ADD_TEST ("reads/test_parse_fq.cpp")

# counting suffixes

#ADD_TEST ("io/test_seq_iter.cpp")
#ADD_TEST ("index/test_suffix.cpp")

# k-mer and terminal extraction

#ADD_TEST ("index/test_kmer_buffer.cpp")
#ADD_TEST ("index/test_terminals.cpp")
#ADD_TEST ("index/test_sorting.cpp")

## flush

#ADD_TEST ("index/test_bit_vectors.cpp")
#ADD_TEST ("index/test_graph.cpp")
#ADD_TEST ("index/test_graph_nav.cpp")
#ADD_TEST ("main/test_genbank.cpp")

## distances

#ADD_TEST ("dist/test_disk_map.cpp")
#ADD_TEST ("dist/test_dist.cpp")

## classify

#ADD_TEST ("dist/test_classify.cpp")

## trees

#ADD_TEST ("tree/test_newick.cpp")

# append executables in parent scope
set (MAKI_EXECUTABLES ${MAKI_EXECUTABLES} PARENT_SCOPE)
