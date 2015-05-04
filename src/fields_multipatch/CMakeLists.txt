ADD_LIBRARY( sll_scalar_field_multipatch STATIC
             sll_scalar_field_2d_multipatch.F90 )

TARGET_LINK_LIBRARIES( sll_scalar_field_multipatch
                       sll_coordinate_transformation_multipatch
                       sll_fields
                       sll_meshes_multipatch )

# Ctest
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_scalar_field_multipatch_2d unit_test_2d.F90 )
  TARGET_LINK_LIBRARIES( test_scalar_field_multipatch_2d sll_scalar_field_multipatch )

  IF(PYTHON3_FOUND)
    ADD_TEST( NAME scalar_field_multipatch_2d COMMAND test_scalar_field_multipatch_2d )
    SET_TESTS_PROPERTIES( scalar_field_multipatch_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(PYTHON3_FOUND)
  
  SET( INPUT_FILES_DIR ${CMAKE_SOURCE_DIR}/CAID_coordinate_transformations/square_4p_n10/ )
  FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ENDIF(BUILD_TESTING)
