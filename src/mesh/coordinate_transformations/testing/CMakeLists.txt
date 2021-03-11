set(INPUT_FILES_DIR
    ${CMAKE_SOURCE_DIR}/add_ons/multipatch/testing/CAID_coordinate_transformations/mesh_collela/
)
file(GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml")
file(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_coordinate_transformations_2d
               test_coordinate_transformations_2d.F90)
target_link_libraries(test_coordinate_transformations_2d
                      sll_coordinate_transformations)
add_test(NAME coordinate_transformations
         COMMAND test_coordinate_transformations_2d)
set_tests_properties(coordinate_transformations
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

# -------------------------------------------------------------------------------
# Unit tests for singular mapped disk-like domains
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# OBJECT libraries containing analytical profiles and error bounds
# -------------------------------------------------------------------------------
add_library(
  analytical_profiles_ct OBJECT
  ../../../splines/splines_basic/testing/m_analytical_profiles_1d.F90
  ../../../splines/splines_basic/testing/m_analytical_profiles_1d_base.F90
  ../../../splines/splines_basic/testing/m_analytical_profiles_1d_cos.F90
  ../../../splines/splines_basic/testing/m_analytical_profiles_1d_poly.F90
  ../../../splines/splines_basic/testing/m_analytical_profiles_2d.F90
  ../../../splines/splines_basic/testing/m_analytical_profiles_2d_base.F90
  ../../../splines/splines_basic/testing/m_analytical_profiles_2d_cos_cos.F90
  ../../../splines/splines_basic/testing/m_analytical_profiles_2d_poly.F90)
add_dependencies(analytical_profiles_ct sll_working_precision sll_constants)

# Test singular disk-like mappings
add_executable(test_singular_mapping test_singular_mapping.F90
                                     $<TARGET_OBJECTS:analytical_profiles_ct>)
target_link_libraries(test_singular_mapping sll_coordinate_transformations
                      sll_splines)

# Test advection on singular disk-like mappings using pseudo-Cartesian
# coordinates
set(TEST_NAME test_singular_mapping_advection)
file(GLOB INPUT_FILE ${TEST_NAME}.nml)
file(COPY ${INPUT_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_executable(${TEST_NAME} ${TEST_NAME}.F90)
target_link_libraries(${TEST_NAME} sll_coordinate_transformations sll_splines)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
                                             TIMEOUT 60)

set(TEST_NAME test_jacobian_2d_pseudo_cartesian)
add_executable(${TEST_NAME} ${TEST_NAME}.F90)
target_link_libraries(${TEST_NAME} sll_coordinate_transformations sll_splines)

# add_executable( test_greville_points test_greville_points.F90 )
# target_link_libraries( test_greville_points sll_splines )
