file(GLOB IO_SRC "sll_*.F90")
add_library(sll_low_level_file_io STATIC ${IO_SRC})
target_link_libraries(sll_low_level_file_io ${HDF5_LIBRARIES})

add_executable(test_io unit_test.F90)
target_link_libraries(test_io sll_low_level_file_io sll_memory sll_assertion)
