# Save .mod files locally (this does not pollute the main "modules" directory)
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# OBJECT libraries containing analytical profiles and error bounds
#-------------------------------------------------------------------------------
add_library( analytical_profiles OBJECT
  m_analytical_profiles_1d.F90
  m_analytical_profiles_1d_base.F90
  m_analytical_profiles_1d_cos.F90
  m_analytical_profiles_1d_poly.F90
  m_analytical_profiles_2d.F90
  m_analytical_profiles_2d_base.F90
  m_analytical_profiles_2d_cos_cos.F90
  m_analytical_profiles_2d_poly.F90 )
add_dependencies( analytical_profiles sll_working_precision sll_constants )

add_library( splines_error_bounds OBJECT m_splines_error_bounds.F90 )
add_dependencies( splines_error_bounds analytical_profiles )

#-------------------------------------------------------------------------------
# New test for new interface of 1D splines
#-------------------------------------------------------------------------------
add_executable( test_spline_1d test_spline_1d.F90 m_test_spline_1d.F90
  $<TARGET_OBJECTS:analytical_profiles>
  $<TARGET_OBJECTS:splines_error_bounds> )
target_link_libraries( test_spline_1d sll_splines )

add_test( NAME spline_1d_uniform                     COMMAND test_spline_1d -u )
add_test( NAME spline_1d_non_uniform_on_uniform_grid COMMAND test_spline_1d -n )
add_test( NAME spline_1d_non_uniform_on_general_grid COMMAND test_spline_1d -n 0.1 )
set_tests_properties(
  spline_1d_uniform
  spline_1d_non_uniform_on_uniform_grid
  spline_1d_non_uniform_on_general_grid
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

#-------------------------------------------------------------------------------
# New tests for 2D tensor-product spline using new interface
#-------------------------------------------------------------------------------
add_executable( test_spline_2d test_spline_2d.F90 m_test_spline_2d.F90
  $<TARGET_OBJECTS:analytical_profiles>
  $<TARGET_OBJECTS:splines_error_bounds> )
target_link_libraries( test_spline_2d sll_splines )

add_test( NAME spline_2d_uniform                     COMMAND test_spline_2d -u )
add_test( NAME spline_2d_non_uniform_on_uniform_grid COMMAND test_spline_2d -n )
add_test( NAME spline_2d_non_uniform_on_general_grid COMMAND test_spline_2d -n 0.1 )
set_tests_properties(
  spline_2d_uniform
  spline_2d_non_uniform_on_uniform_grid
  spline_2d_non_uniform_on_general_grid
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

#-------------------------------------------------------------------------------
# New tests for 2D polar spline (experimental)
#-------------------------------------------------------------------------------
add_executable( test_polar_spline_2d test_polar_spline_2d.F90
  $<TARGET_OBJECTS:analytical_profiles> )
target_link_libraries( test_polar_spline_2d sll_splines )

#-------------------------------------------------------------------------------

# Cubic splines ............................
ADD_EXECUTABLE( test_cubic_splines test_cubic_splines.F90  )
TARGET_LINK_LIBRARIES( test_cubic_splines sll_splines )

# Ctests
ADD_TEST( NAME cubic_splines COMMAND test_cubic_splines )
SET_TESTS_PROPERTIES( cubic_splines PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................

# Quintic splines ............................
ADD_EXECUTABLE( test_quintic_splines test_quintic_splines.F90 )
TARGET_LINK_LIBRARIES( test_quintic_splines sll_splines )

# Ctests
ADD_TEST( NAME quintic_splines COMMAND test_quintic_splines )
SET_TESTS_PROPERTIES( quintic_splines PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................

# Cubic non uniform splines ............................
ADD_EXECUTABLE( test_cubic_splines_nonuniform test_cubic_splines_nonuniform.F90)
TARGET_LINK_LIBRARIES( test_cubic_splines_nonuniform sll_splines sll_linear_solvers )

# Ctests
ADD_TEST( NAME cubic_splines_nonuniform COMMAND test_cubic_splines_nonuniform )
SET_TESTS_PROPERTIES( cubic_splines_nonuniform PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................

# Basis splines ............................
ADD_EXECUTABLE( test_bsplines test_bsplines.F90)
TARGET_LINK_LIBRARIES( test_bsplines sll_splines )

# Ctests
ADD_TEST( NAME bsplines  COMMAND test_bsplines )
SET_TESTS_PROPERTIES( bsplines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................

# Box-splines ............................
ADD_EXECUTABLE( test_box_splines_derivatives test_box_splines_derivatives.F90 )
TARGET_LINK_LIBRARIES( test_box_splines_derivatives sll_splines )


# Splines_pp...........................................
ADD_EXECUTABLE( test_splines_pp_1d test_splines_pp_1d.F90 )
TARGET_LINK_LIBRARIES( test_splines_pp_1d sll_splines sll_interpolators )
#...........................................
ADD_EXECUTABLE( test_splines_pp_2d test_splines_pp_2d.F90 )
TARGET_LINK_LIBRARIES( test_splines_pp_2d sll_splines sll_interpolators )
#...........................................
ADD_EXECUTABLE( test_splines_pp_3d test_splines_pp_3d.F90 )
TARGET_LINK_LIBRARIES( test_splines_pp_3d sll_splines sll_interpolators )

# Ctests
ADD_TEST( NAME splines_pp_1d COMMAND test_splines_pp_1d)
ADD_TEST( NAME splines_pp_2d COMMAND test_splines_pp_2d) 
ADD_TEST( NAME splines_pp_3d COMMAND test_splines_pp_3d) 
SET_TESTS_PROPERTIES( splines_pp_1d PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( splines_pp_2d PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( splines_pp_3d PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................
