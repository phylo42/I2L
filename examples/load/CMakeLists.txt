cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Create an executable target and set properties
add_executable(i2l_example_load EXCLUDE_FROM_ALL "")

target_sources(i2l_example_load
    PRIVATE
        src/main.cpp
)

find_package(Boost REQUIRED COMPONENTS filesystem)

target_link_libraries(i2l_example_load
    PRIVATE
        i2l::dna
        Boost::filesystem
)

# Turn on the warnings and treat them as errors
target_compile_options(i2l_example_load
    PRIVATE
        -Wall -Wextra -Wpedantic -Werror
)

target_compile_features(i2l_example_load
    PRIVATE
        cxx_std_17
)