set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_library(sll_parallel_array_initializers STATIC 
  sll_parallel_array_initializer_module.F90
  sll_common_array_initializers_module.F90 
  )

add_dependencies(sll_parallel_array_initializers
  sll_remap
  sll_coordinate_transformations
  sll_logical_meshes
  sll_assert
  sll_working_precision
  )

add_executable(test_parallel_array_initializer unit_test_4d.F90)
target_link_libraries(test_parallel_array_initializer
	sll_remap
	sll_collective
	sll_coordinate_transformations
	sll_file_io
	sll_logical_meshes
	sll_memory
	sll_assertion
	sll_working_precision
        ${HDF5_LIBRARIES}
	)

# Just in case that the unit test needs diagnostics
#IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

#   ADD_DEPENDENCIES(sll_remap sll_low_level_file_io_parallel)
#   TARGET_LINK_LIBRARIES(sll_remap ${HDF5_LIBRARIES})
#   ADD_EXECUTABLE(test_layout_output unit_test_output.F90)
#   TARGET_LINK_LIBRARIES(test_layout_output sll_remap
#                                            sll_low_level_file_io_parallel
#                                            sll_collective
#	                                    sll_utilities
#	                                    sll_constants
#	                                    sll_memory
#	                                    sll_assertion
#	                                    sll_working_precision)




