if(HDF5_ENABLED
   AND HDF5_PARALLEL_ENABLED
   AND HDF5_IS_PARALLEL)
  add_executable(
    test_parallel_array_initializer test_parallel_array_initializer.F90
                                    m_parallel_array_output.F90)

  target_link_libraries(test_parallel_array_initializer
                        sll_parallel_array_utilities)

  set(ARGS " ")
  add_mpi_test(parallel_array_initializers test_parallel_array_initializer
               ${PROCS} ${ARGS})
  set_tests_properties(parallel_array_initializers
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

endif()
