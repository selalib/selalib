ADD_EXECUTABLE( test_particle_visualization_interface test_particle_visualization_interface.F90)
TARGET_LINK_LIBRARIES( test_particle_visualization_interface
  sll_pic_visualization
  sll_pic_sampling
  sll_pic_weights
  sll_xdmf_io
  )

SET( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/distribution_for_the_test.nml )
ADD_TEST( NAME particle_visualization_interface COMMAND test_particle_visualization_interface ${NML_FILE})
SET_TESTS_PROPERTIES(particle_visualization_interface PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
