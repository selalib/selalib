add_executable(test_particle_group_1d1v test_particle_group_1d1v.F90)
target_link_libraries(test_particle_group_1d1v sll_particle_groups)

add_executable(test_particle_group_1d2v test_particle_group_1d2v.F90)
target_link_libraries(test_particle_group_1d2v sll_particle_groups)

add_executable(test_particle_group_2d2v test_particle_group_2d2v.F90)
target_link_libraries(test_particle_group_2d2v sll_particle_groups)

add_executable(test_particle_group_3d1v test_particle_group_3d1v.F90)
target_link_libraries(test_particle_group_3d1v sll_particle_groups)

add_test(NAME particle_group_1d1v COMMAND test_particle_group_1d1v)
set_tests_properties(particle_group_1d1v PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

add_test(NAME particle_group_1d2v COMMAND test_particle_group_1d2v)
set_tests_properties(particle_group_1d2v PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

add_test(NAME particle_group_2d2v COMMAND test_particle_group_2d2v)
set_tests_properties(particle_group_2d2v PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

add_test(NAME particle_group_3d1v COMMAND test_particle_group_3d1v)
set_tests_properties(particle_group_3d1v PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

if(NOT (FFT_LIB STREQUAL "SLLFFT"))

  add_executable(test_particle_group_2d2v_lbf test_particle_group_2d2v_lbf.F90)
  target_link_libraries(test_particle_group_2d2v_lbf sll_particle_groups)

  add_test(NAME particle_group_2d2v_lbf COMMAND test_particle_group_2d2v_lbf)
  set_tests_properties(particle_group_2d2v_lbf
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

endif()
