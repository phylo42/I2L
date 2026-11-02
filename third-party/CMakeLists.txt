cmake_minimum_required(VERSION 4.1)

execute_process(
        COMMAND patch -p1 -i ${CMAKE_CURRENT_SOURCE_DIR}/hopscotch_cmake_update.patch
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/hopscotch-map
        RESULT_VARIABLE patch_result
)

if(NOT patch_result EQUAL 0)
    message(STATUS "Patch already applied or failed (this is usually normal after the first run)")
endif()

add_subdirectory(hopscotch-map)
