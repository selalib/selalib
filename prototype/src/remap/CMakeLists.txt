set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_library(sll_remap STATIC sll_remap.F90)
add_dependencies(sll_remap sll_collective sll_memory sll_utilities sll_working_precision)
# test remap 3d
add_executable(test_remap unit_test.F90)
target_link_libraries(test_remap
	sll_remap
	sll_utilities
	sll_collective
	sll_memory
	sll_assertion
	sll_working_precision
)

#test remap 2d
add_executable(test_remap_2d unit_test_2d.F90)
target_link_libraries(test_remap_2d
	sll_remap
	sll_collective
	sll_utilities
	sll_memory
	sll_assertion
	sll_working_precision
)

#test remap 4d
add_executable(test_remap_4d unit_test_4d.F90)
target_link_libraries(test_remap_4d
	sll_remap
	sll_collective
	sll_utilities
	sll_memory
	sll_assertion
	sll_working_precision
)



IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)

   ADD_DEPENDENCIES(sll_remap sll_low_level_file_io_parallel)
   TARGET_LINK_LIBRARIES(sll_remap ${HDF5_LIBRARIES})
   ADD_EXECUTABLE(test_layout_output unit_test_output.F90)
   TARGET_LINK_LIBRARIES(test_layout_output sll_remap
                                            sll_low_level_file_io_parallel
                                            sll_collective
	                                    sll_utilities
	                                    sll_constants
	                                    sll_memory
	                                    sll_assertion
	                                    sll_working_precision)
ENDIF()
