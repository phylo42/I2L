cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

add_subdirectory(amino_acids)
add_subdirectory(basic)
add_subdirectory(load)
add_subdirectory(search)
add_subdirectory(stats)
add_subdirectory(kmers)
add_subdirectory(tree)

add_custom_target(examples
        DEPENDS
            i2l_example_aa
            i2l_example_basic
            i2l_example_load
            i2l_example_search
            i2l_example_stats
            i2l_example_kmers
            i2l_example_tree)