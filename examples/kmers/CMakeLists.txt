cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Create an executable target and set properties
add_executable(i2l_example_kmers EXCLUDE_FROM_ALL "")

target_sources(i2l_example_kmers
        PRIVATE
            src/main.cpp
        )

target_link_libraries(i2l_example_kmers
    PRIVATE
        i2l::dna
)

# Turn on the warnings and treat them as errors
target_compile_options(i2l_example_kmers
    PRIVATE
        -Wall -Wextra -Wpedantic -Werror
)

target_compile_features(i2l_example_kmers
    PRIVATE
        cxx_std_17
)