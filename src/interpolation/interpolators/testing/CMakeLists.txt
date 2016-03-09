SLL_TEST( cubic_spline_interpolator_1d sll_interpolators )
SLL_TEST( cubic_spline_interpolator_2d sll_interpolators )
SLL_TEST( arbitrary_degree_spline_interpolator_1d sll_interpolators sll_file_io )
SLL_TEST( arbitrary_degree_spline_interpolator_2d sll_interpolators sll_file_io )
SLL_TEST( periodic_interpolator_1d sll_interpolators )
SLL_TEST( arbitrary_degree_spline_1d_periodic  sll_interpolators )
SLL_TEST( arbitrary_degree_spline_1d_dirichlet sll_interpolators )
SLL_TEST( arbitrary_degree_spline_1d_neumann sll_interpolators )
SLL_TEST( quintic_spline_interpolator_1d sll_interpolators )
#SLL_TEST( bspline_interpolator_1d sll_interpolators )
#SLL_TEST( bspline_interpolator_2d sll_interpolators )

SET_TESTS_PROPERTIES( cubic_spline_interpolator_1d 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
SET_TESTS_PROPERTIES( cubic_spline_interpolator_2d 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
SET_TESTS_PROPERTIES( arbitrary_degree_spline_interpolator_1d 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
SET_TESTS_PROPERTIES( arbitrary_degree_spline_interpolator_2d 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )
SET_TESTS_PROPERTIES( periodic_interpolator_1d 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
SET_TESTS_PROPERTIES( arbitrary_degree_spline_1d_periodic
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
SET_TESTS_PROPERTIES( arbitrary_degree_spline_1d_dirichlet
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
SET_TESTS_PROPERTIES( arbitrary_degree_spline_1d_neumann 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
SET_TESTS_PROPERTIES( quintic_spline_interpolator_1d 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 5 )
