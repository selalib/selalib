set(ARGS " ")

if(PROCESSOR_COUNT GREATER 1)

  add_executable(test_remap_2d test_remap_2d.F90)
  target_link_libraries(test_remap_2d sll_remap)
  add_mpi_test(remap_2d test_remap_2d 2 ${ARGS})
  set_tests_properties(remap_2d PROPERTIES TIMEOUT 20 PASS_REGULAR_EXPRESSION
                                           "PASSED")

  add_executable(test_remap_3d test_remap_3d.F90)
  target_link_libraries(test_remap_3d sll_remap)
  add_mpi_test(remap_3d test_remap_3d 2 ${ARGS})
  set_tests_properties(remap_3d PROPERTIES TIMEOUT 30 PASS_REGULAR_EXPRESSION
                                           "PASSED")

  add_executable(test_remap_4d test_remap_4d.F90)
  target_link_libraries(test_remap_4d sll_remap)
  add_mpi_test(remap_4d test_remap_4d 2 ${ARGS})
  set_tests_properties(remap_4d PROPERTIES TIMEOUT 60 PASS_REGULAR_EXPRESSION
                                           "PASSED")

  add_executable(test_remap_5d test_remap_5d.F90)
  target_link_libraries(test_remap_5d sll_remap)
  add_mpi_test(remap_5d test_remap_5d 2 ${ARGS})
  set_tests_properties(remap_5d PROPERTIES TIMEOUT 120 PASS_REGULAR_EXPRESSION
                                           "PASSED")

  add_executable(test_remap_6d test_remap_6d.F90)
  target_link_libraries(test_remap_6d sll_remap)
  add_mpi_test(remap_6d test_remap_6d 2 ${ARGS})
  set_tests_properties(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
                                           TIMEOUT 240)

  # ADD_EXECUTABLE( test_remap_6d_transp test_remap_6d_transp.F90 )
  # TARGET_LINK_LIBRARIES( test_remap_6d_transp sll_remap ) ADD_MPI_TEST(
  # remap_6d_transp test_remap_6d_transp ${PROCS} ${ARGS} )
  # SET_TESTS_PROPERTIES( remap_6d_transp PROPERTIES PASS_REGULAR_EXPRESSION
  # "PASSED" TIMEOUT 200 )

endif(PROCESSOR_COUNT GREATER 1)

if(HDF5_ENABLED
   AND HDF5_IS_PARALLEL
   AND HDF5_PARALLEL_ENABLED)
  # Layout output (not in Ctest)
  add_executable(test_layout_output test_layout_output.F90)
  target_link_libraries(test_layout_output sll_remap sll_file_io_parallel)

  # Parallel I/O
  add_executable(test_io_parallel test_io_parallel.F90)
  target_link_libraries(test_io_parallel sll_file_io_parallel sll_remap)
  set(PROCS 4)
  add_mpi_test(io_parallel test_io_parallel ${PROCS} ${ARGS})
  set_tests_properties(io_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
                                              TIMEOUT 20)

endif()
