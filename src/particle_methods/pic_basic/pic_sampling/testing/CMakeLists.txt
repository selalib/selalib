add_executable(test_particle_sampling test_particle_sampling.F90)
target_link_libraries(test_particle_sampling sll_pic_sampling sll_pic_weights)

add_test(NAME particle_sampling COMMAND test_particle_sampling)
set_tests_properties(particle_sampling PROPERTIES PASS_REGULAR_EXPRESSION
                                                  "PASSED")

if(NOT (FFT_LIB STREQUAL "SLLFFT"))

  add_executable(test_particle_sampling_interface
                 test_particle_sampling_interface.F90)
  target_link_libraries(test_particle_sampling_interface sll_pic_sampling
                        sll_particle_groups sll_pic_weights)

  set(NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/sampling_interface_test.nml)
  add_test(NAME particle_sampling_interface
           COMMAND test_particle_sampling_interface ${NML_FILE})
  set_tests_properties(particle_sampling_interface
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif()
