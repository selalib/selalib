set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_library(sll_remap STATIC sll_remap.F90)
add_dependencies(sll_remap 
  sll_collective 
  sll_electric_field_accumulators 
  sll_memory 
  sll_utilities 
  sll_working_precision
  )

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