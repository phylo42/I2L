cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

enable_testing()

add_subdirectory(catch2 EXCLUDE_FROM_ALL)

add_executable(i2l_tests EXCLUDE_FROM_ALL src/main.cpp)
find_package(Boost REQUIRED COMPONENTS filesystem iostreams)

target_link_libraries(i2l_tests
        PRIVATE
            i2l::dna
            Catch2::Catch2WithMain
            Boost::filesystem
            Boost::iostreams
        )

include(CTest)
include(Catch)
catch_discover_tests(i2l_tests)