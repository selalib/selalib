ADD_LIBRARY( sll_bsl_lt_pic STATIC
  bsl_lt_pic_4d_particle.F90
  bsl_lt_pic_4d_group.F90
  bsl_lt_pic_4d_utilities.F90
  )

TARGET_LINK_LIBRARIES( sll_bsl_lt_pic
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes
  # sll_pic_interface
  sll_remapped_pic
  sll_remapped_pic_utilities
  sll_pic_accumulators
  sll_timer # needed?
  sll_file_io # needed?
  )


ADD_EXECUTABLE(test_bsl_lt_pic_4d_remap unit_test_bsl_lt_pic_4d.F90)

TARGET_LINK_LIBRARIES(test_bsl_lt_pic_4d_remap
  sll_timer
  sll_constants
  sll_working_precision
  sll_assert
  sll_remapped_pic
  sll_bsl_lt_pic
  )


# Ctest
IF(BUILD_TESTING)
  # TODO
  # ADD_EXECUTABLE( test_bsl_lt_pic unit_test.F90 )
  # TARGET_LINK_LIBRARIES( test_bsl_lt_pic sll_bsl_lt_pic_pic )
  # ADD_TEST( NAME bsl_lt_pic COMMAND test_bsl_lt_pic )
  # SET_TESTS_PROPERTIES( bsl_lt_pic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)
