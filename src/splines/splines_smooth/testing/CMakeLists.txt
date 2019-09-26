# Save .mod files locally (this does not pollute the main "modules" directory)
#SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Unit tests for singular mapped disk-like domains
#-------------------------------------------------------------------------------

# Test singular disk-like mappings
add_executable( test_singular_mapping test_singular_mapping.F90
  $<TARGET_OBJECTS:analytical_profiles> )
target_link_libraries( test_singular_mapping sll_splines_smooth )

# Test advection on singular disk-like mappings using pseudo-Cartesian coordinates
set( TEST_NAME test_singular_mapping_advection )
file( GLOB INPUT_FILE ${TEST_NAME}.nml )
file( COPY ${INPUT_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
add_executable( ${TEST_NAME} ${TEST_NAME}.F90 )
target_link_libraries( ${TEST_NAME} sll_splines_smooth )
add_test( NAME ${TEST_NAME} COMMAND ${TEST_NAME} )
set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 60 )

set( TEST_NAME test_jacobian_2d_pseudo_cartesian )
add_executable( ${TEST_NAME} ${TEST_NAME}.F90 )
target_link_libraries( ${TEST_NAME} sll_splines_smooth )

#add_executable( test_greville_points test_greville_points.F90 )
#target_link_libraries( test_greville_points sll_splines )
