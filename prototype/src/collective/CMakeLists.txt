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

# we can call ${MPIEXEC_MAX_NUMPROCS} to have the number of procs
SET(PROCS 2)
SET(ARGS "")
ADD_TEST(NAME collective
         COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_collective
	${MPIEXEC_POSTFLAGS} ${ARGS})
SET_TESTS_PROPERTIES(collective PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS")
