# Save .mod files locally (this does not pollute the main "modules" directory)
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# New tests for 2D B-splines
#-------------------------------------------------------------------------------

add_executable( test_bsplines_2d_new
  test_bsplines_2d_new.F90
  m_test_bsplines_2d.F90
  m_analytical_profiles_2d.F90
  m_analytical_profiles_2d_base.F90
  m_analytical_profiles_2d_cos_cos.F90
  m_analytical_profiles_2d_poly.F90
  )
target_link_libraries( test_bsplines_2d_new sll_splines )

# WARNING: test name is not consistent with executable (missing 's') -> TODO
add_test( NAME bspline_2d_new COMMAND test_bsplines_2d_new )
set_tests_properties( bspline_2d_new PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

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

# B-splines ............................
ADD_EXECUTABLE( test_bspline_1d test_bsplines_1d.F90 )
TARGET_LINK_LIBRARIES( test_bspline_1d sll_splines )
ADD_EXECUTABLE( test_bspline_2d test_bsplines_2d.F90 )
TARGET_LINK_LIBRARIES( test_bspline_2d sll_splines )

# Ctests
ADD_TEST( NAME bspline_1d COMMAND test_bspline_1d )
ADD_TEST( NAME bspline_2d COMMAND test_bspline_2d )
SET_TESTS_PROPERTIES( bspline_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( bspline_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
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
