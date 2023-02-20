cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

add_subdirectory(hopscotch-map)
set(HASH_MAP_IMPLEMENTATION tsl::hopscotch_map)
