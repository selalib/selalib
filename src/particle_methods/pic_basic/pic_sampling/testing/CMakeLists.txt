ADD_EXECUTABLE( test_particle_sampling test_particle_sampling.F90)
TARGET_LINK_LIBRARIES( test_particle_sampling sll_pic_sampling sll_pic_weights)

ADD_TEST( NAME particle_sampling COMMAND test_particle_sampling)
SET_TESTS_PROPERTIES(particle_sampling PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

if (NOT (FFT_LIB STREQUAL "SLLFFT"))
 
  ADD_EXECUTABLE( test_particle_sampling_interface test_particle_sampling_interface.F90)
  TARGET_LINK_LIBRARIES( test_particle_sampling_interface sll_pic_sampling sll_particle_groups sll_pic_weights)

  SET( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/sampling_interface_test.nml )
  ADD_TEST( NAME particle_sampling_interface COMMAND test_particle_sampling_interface ${NML_FILE})
  SET_TESTS_PROPERTIES(particle_sampling_interface PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif()
