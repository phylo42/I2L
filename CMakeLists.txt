cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(I2L LANGUAGES CXX VERSION 0.4.0)

find_package(Boost REQUIRED COMPONENTS serialization iostreams system filesystem)
find_package(ZLIB)

add_subdirectory(third-party)
add_subdirectory(tests)
add_subdirectory(examples)

set(SOURCES
        src/i2l/fasta.cpp include/i2l/fasta.h
        src/i2l/file_io.cpp include/i2l/file_io.h
        include/i2l/hash_map.h
        src/i2l/kmer_iterator.cpp include/i2l/kmer_iterator.h
        src/i2l/newick.cpp include/i2l/newick.h
        include/i2l/optional.h
        src/i2l/phylo_kmer.cpp include/i2l/phylo_kmer.h
        include/i2l/phylo_kmer_db.h
        src/i2l/phylo_node.cpp include/i2l/phylo_node.h
        src/i2l/phylo_tree.cpp include/i2l/phylo_tree.h
        src/i2l/seq.cpp include/i2l/seq.h
        src/i2l/seq_record.cpp include/i2l/seq_record.h
        include/i2l/serialization.h src/i2l/serialization.cpp
        include/i2l/version.h)


set(LINK_LIBRARIES
        Boost::serialization Boost::iostreams Boost::system Boost::filesystem
        ZLIB::ZLIB
        tsl::hopscotch_map)

set(INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
set(COMPILE_OPTIONS -Wall -Wextra -Wpedantic)

# explicitly require c++17 standard
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
set(COMPILE_FEATURES cxx_std_17)

##########################################################################################
# Proteins
add_library(i2l_aa STATIC "")
add_library(i2l::aa ALIAS i2l_aa)

target_sources(i2l_aa PRIVATE ${SOURCES})
target_link_libraries(i2l_aa PUBLIC ${LINK_LIBRARIES})
target_include_directories(i2l_aa PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(i2l_aa PRIVATE ${COMPILE_OPTIONS})
target_compile_features(i2l_aa PUBLIC ${COMPILE_FEATURES})
target_compile_definitions(i2l_aa
        PUBLIC
        SEQ_TYPE_AA)

###########################################################################################
## Proteins, positions included
add_library(i2l_aa_pos STATIC "")
add_library(i2l::aa_pos ALIAS i2l_aa_pos)

target_sources(i2l_aa_pos PRIVATE ${SOURCES})
target_link_libraries(i2l_aa_pos PUBLIC ${LINK_LIBRARIES})
target_include_directories(i2l_aa_pos PUBLIC ${INCLUDE_DIRECTORIES})
set_property(TARGET i2l_aa_pos PROPERTY CXX_STANDARD 17)
target_compile_options(i2l_aa_pos PRIVATE ${COMPILE_OPTIONS})
target_compile_features(i2l_aa_pos PUBLIC ${COMPILE_FEATURES})
target_compile_definitions(i2l_aa_pos
        PUBLIC
            SEQ_TYPE_AA
            KEEP_POSITIONS)

##########################################################################################
# DNA
# Create a library target and set properties
add_library(i2l_dna STATIC "")
add_library(i2l::dna ALIAS i2l_dna)

target_sources(i2l_dna PRIVATE ${SOURCES})
target_link_libraries(i2l_dna PUBLIC ${LINK_LIBRARIES})
target_include_directories(i2l_dna PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(i2l_dna PRIVATE ${COMPILE_OPTIONS})
target_compile_features(i2l_dna PUBLIC ${COMPILE_FEATURES})
set_property(TARGET i2l_dna PROPERTY CXX_STANDARD 17)
target_compile_definitions(i2l_dna
        PUBLIC
            SEQ_TYPE_DNA)

######################################################################################################
# Library installation
include(GNUInstallDirs)
#include(CMakePackageConfigHelpers)

# Exporting the core target
# 'make install' or `cmake --install` to the correct locations
# (provided by GNUInstallDirs)
install(TARGETS i2l_dna
        EXPORT i2lDnaTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS i2l_aa
        EXPORT i2lAaTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#export(TARGETS i2l_dna i2l_aa ${HASH_MAP_EXPORT}
export(TARGETS i2l_dna hopscotch_map #${HASH_MAP_EXPORT}
       FILE ${CMAKE_BINARY_DIR}/i2lConfig.cmake)

export(PACKAGE i2l)
