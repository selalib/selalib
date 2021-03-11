add_executable(test_point_to_point_comms_1d test_point_to_point_comms_1d.F90)
target_link_libraries(
  test_point_to_point_comms_1d sll_point_to_point_communications sll_collective
  sll_timer)

add_executable(test_point_to_point_comms_2d test_point_to_point_comms_2d.F90)
target_link_libraries(
  test_point_to_point_comms_2d sll_point_to_point_communications sll_collective
  sll_timer)

set(ARGS " ")
add_mpi_test(point_to_point_comms_1d test_point_to_point_comms_1d ${PROCS}
             ${ARGS})
add_mpi_test(point_to_point_comms_2d test_point_to_point_comms_2d "4" ${ARGS})
set_tests_properties(point_to_point_comms_1d PROPERTIES PASS_REGULAR_EXPRESSION
                                                        "PASSED")
set_tests_properties(point_to_point_comms_2d PROPERTIES PASS_REGULAR_EXPRESSION
                                                        "PASSED")
