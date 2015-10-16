SET( INPUT_FILES_DIR ${CMAKE_SOURCE_DIR}/add_ons/multipatch/testing/CAID_coordinate_transformations/mesh_collela/ )
FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ADD_EXECUTABLE( test_coordinate_transformations_2d test_coordinate_transformations_2d.F90)
TARGET_LINK_LIBRARIES( test_coordinate_transformations_2d sll_coordinate_transformations )
ADD_TEST( NAME coordinate_transformations COMMAND test_coordinate_transformations_2d )
SET_TESTS_PROPERTIES( coordinate_transformations PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
