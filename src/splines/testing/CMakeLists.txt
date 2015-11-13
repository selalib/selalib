# Cubic splines ............................
ADD_EXECUTABLE( test_cubic_splines test_cubic_splines.F90 test_func_module.F90
                util_constants.F90 test_processes_module.F90 )
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

# Deboor B-splines ............................
ADD_EXECUTABLE( test_bsplines_1d test_bsplines_1d.F90 )
TARGET_LINK_LIBRARIES( test_bsplines_1d sll_splines )

ADD_EXECUTABLE( test_bsplines_2d test_bsplines_2d.F90 )
TARGET_LINK_LIBRARIES( test_bsplines_2d sll_splines )

# Ctests
ADD_TEST( NAME bsplines_1d COMMAND test_bsplines_1d )
ADD_TEST( NAME bsplines_2d COMMAND test_bsplines_2d )
SET_TESTS_PROPERTIES( bsplines_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( bsplines_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................

# Cubic non uniform splines ............................
ADD_EXECUTABLE( test_cubic_splines_nonuniform test_cubic_splines_nonuniform.F90)
TARGET_LINK_LIBRARIES( test_cubic_splines_nonuniform sll_splines sll_linear_solvers )

# Ctests
ADD_TEST( NAME cubic_splines_nonuniform COMMAND test_cubic_splines_nonuniform )
SET_TESTS_PROPERTIES( cubic_splines_nonuniform PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................

# Arbitrary degree splines ............................
ADD_EXECUTABLE( test_arbitrary_degree_splines test_arbitrary_degree_splines.F90)
TARGET_LINK_LIBRARIES( test_arbitrary_degree_splines sll_splines )

# Ctests
ADD_TEST( NAME arbitrary_degree_splines  COMMAND test_arbitrary_degree_splines )
SET_TESTS_PROPERTIES( arbitrary_degree_splines PROPERTIES
                      PASS_REGULAR_EXPRESSION "PASSED" )
#...........................................

# Box-splines ............................
ADD_EXECUTABLE( test_box_splines_derivatives test_box_splines_derivatives.F90 )
TARGET_LINK_LIBRARIES( test_box_splines_derivatives sll_splines )

# Ctests
#...........................................
FIND_PACKAGE(GSL)
IF(GSL_FOUND)
  ADD_EXECUTABLE( test_gsl_bsplines test_gsl_bsplines.c )
  TARGET_LINK_LIBRARIES( test_gsl_bsplines  ${GSL_LIB})
ENDIF(GSL_FOUND)
