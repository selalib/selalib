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

SET(PROCS 8)
SET(ARGS "")
ADD_TEST(NAME remap
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap
	${MPIEXEC_POSTFLAGS} ${ARGS})
SET_TESTS_PROPERTIES(remap PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

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

ADD_TEST(NAME remap_2d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_2d
	${MPIEXEC_POSTFLAGS} ${ARGS})
SET_TESTS_PROPERTIES(remap_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

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

SET(PROCS 16)
ADD_TEST(NAME remap_4d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_4d
	${MPIEXEC_POSTFLAGS} ${ARGS})
SET_TESTS_PROPERTIES(remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

#test remap 6d
add_executable(test_remap_6d unit_test_6d.F90)
target_link_libraries(test_remap_6d
	sll_remap
	sll_collective
	sll_utilities
	sll_memory
	sll_assertion
	sll_working_precision
)

SET(PROCS 16)
ADD_TEST(NAME remap_6d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_6d
	${MPIEXEC_POSTFLAGS} ${ARGS})
SET_TESTS_PROPERTIES(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

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



