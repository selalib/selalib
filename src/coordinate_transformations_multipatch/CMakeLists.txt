ADD_LIBRARY( sll_coordinate_transformation_multipatch STATIC
  sll_coordinate_transformation_multipatch.F90 )

TARGET_LINK_LIBRARIES( sll_coordinate_transformation_multipatch
  sll_coordinate_transformations
  sll_meshes_multipatch
  sll_memory )

# Ctest
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_coordinate_transformation_multipatch_2d unit_test_2d.F90 )
  TARGET_LINK_LIBRARIES(test_coordinate_transformation_multipatch_2d
     sll_coordinate_transformation_multipatch
     sll_file_io
     sll_memory)

  IF(PYTHON3_FOUND)
    ADD_TEST( NAME coordinate_transformation_multipatch_2d
    COMMAND test_coordinate_transformation_multipatch_2d )
    SET_TESTS_PROPERTIES( coordinate_transformation_multipatch_2d 
    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(PYTHON3_FOUND)

  SET( INPUT_FILES_DIR ${CMAKE_SOURCE_DIR}/CAID_coordinate_transformations/circle_mp5_pts12/ )
  FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ENDIF(BUILD_TESTING)
