ADD_EXECUTABLE( test_particle_group_1d1v test_particle_group_1d1v.F90 )
TARGET_LINK_LIBRARIES( test_particle_group_1d1v sll_particle_groups )

ADD_EXECUTABLE( test_particle_group_1d2v test_particle_group_1d2v.F90 )
TARGET_LINK_LIBRARIES( test_particle_group_1d2v sll_particle_groups )

ADD_EXECUTABLE( test_particle_group_2d2v test_particle_group_2d2v.F90 )
TARGET_LINK_LIBRARIES( test_particle_group_2d2v sll_particle_groups )

ADD_EXECUTABLE( test_particle_group_2d3v test_particle_group_2d3v.F90 )
TARGET_LINK_LIBRARIES( test_particle_group_2d3v sll_particle_groups )

ADD_EXECUTABLE( test_particle_group_3d1v test_particle_group_3d1v.F90 )
TARGET_LINK_LIBRARIES( test_particle_group_3d1v sll_particle_groups )

ADD_TEST(NAME particle_group_1d1v COMMAND test_particle_group_1d1v)
SET_TESTS_PROPERTIES(particle_group_1d1v PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME particle_group_1d2v COMMAND test_particle_group_1d2v)
SET_TESTS_PROPERTIES(particle_group_1d2v PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME particle_group_2d2v COMMAND test_particle_group_2d2v)
SET_TESTS_PROPERTIES(particle_group_2d2v PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME particle_group_2d3v COMMAND test_particle_group_2d3v)
SET_TESTS_PROPERTIES(particle_group_2d3v PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME particle_group_3d1v COMMAND test_particle_group_3d1v)
SET_TESTS_PROPERTIES(particle_group_3d1v PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

if (NOT (FFT_LIB STREQUAL "SLLFFT"))
  
  ADD_EXECUTABLE( test_particle_group_2d2v_lbf test_particle_group_2d2v_lbf.F90 )
  TARGET_LINK_LIBRARIES( test_particle_group_2d2v_lbf sll_particle_groups )

  ADD_TEST(NAME particle_group_2d2v_lbf COMMAND test_particle_group_2d2v_lbf)
  SET_TESTS_PROPERTIES(particle_group_2d2v_lbf PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

endif()   
