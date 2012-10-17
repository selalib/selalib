set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_library(sll_collective STATIC sll_collective.F90)
add_dependencies(sll_collective
        sll_utilities
	sll_memory
	sll_assertion
        sll_timer
	sll_working_precision
)

add_executable(test_collective unit_test.F90)
target_link_libraries(test_collective
	sll_collective
	sll_timer
	sll_memory
	sll_assertion
	sll_working_precision
)
