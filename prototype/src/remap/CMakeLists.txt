set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_library(sll_remap STATIC sll_remap.F90)
add_executable(test_remap unit_test.F90)
target_link_libraries(test_remap
	sll_remap
	sll_utilities
	sll_collective
	sll_memory
	sll_assertion
)
