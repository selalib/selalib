SET( INTERPOLATORS_SRC_FILES 
        sll_interpolators_1d_base.F90
        sll_interpolators_2d_base.F90 
        sll_cubic_spline_interpolator_1d.F90
        sll_cubic_spline_interpolator_1d_nonuniform.F90
        sll_quintic_spline_interpolator_1d.F90
        sll_quintic_spline_interpolator_1d_nonuniform.F90
        sll_odd_degree_spline_interpolator_1d.F90
        sll_odd_degree_spline_interpolator_1d_nonuniform.F90
        sll_lagrange_interpolator_1d.F90 
        sll_periodic_interpolator_1d.F90
        sll_cubic_spline_interpolator_2d.F90
        sll_hermite_interpolator_2d.F90
	  sll_arbitrary_degree_spline_interpolator_1d.F90
        sll_arbitrary_degree_spline_interpolator_2d.F90) 

ADD_LIBRARY( sll_interpolators STATIC ${INTERPOLATORS_SRC_FILES} )

TARGET_LINK_LIBRARIES( sll_interpolators sll_splines
                                         sll_periodic_interp
                                         sll_lagrange_interpolation
                                         sll_hermite_interpolation
                                         sll_fft
                                         sll_utilities 
   		                             deboor_splines 
					           sll_timer
                                         ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(test_interpolators_1d unit_test_1d.F90)

TARGET_LINK_LIBRARIES( test_interpolators_1d sll_interpolators)
ADD_TEST(NAME interpolators_1d COMMAND test_interpolators_1d)
SET_TESTS_PROPERTIES(interpolators_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_arb_deg_spline_interpolators_2d 
               unit_test_arbitrary_degree_2d.F90 )

TARGET_LINK_LIBRARIES( test_arb_deg_spline_interpolators_2d
                       sll_interpolators sll_file_io)

ADD_EXECUTABLE(test_arb_deg_spline_interpolators_1d 
               unit_test_arbitrary_degree_1d.F90 )

TARGET_LINK_LIBRARIES( test_arb_deg_spline_interpolators_1d
                          sll_interpolators sll_file_io)

#Depends on coordinate transformations (not possible)
#ADD_EXECUTABLE(test_interpolators_2d unit_test.F90)
#  
#TARGET_LINK_LIBRARIES( test_interpolators_2d sll_interpolators sll_utilities )
#ADD_TEST(NAME interpolators_2d COMMAND test_interpolators_2d)
#SET_TESTS_PROPERTIES(interpolators_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_quintic_interpolators_1d 
               test_sll_quintic_spline_interpolator_1d.F90)

TARGET_LINK_LIBRARIES( test_quintic_interpolators_1d sll_interpolators)

ADD_EXECUTABLE(test_quintic_interpolators_1d_nonuniform 
               test_sll_quintic_spline_interpolator_1d_nonuniform.F90)

TARGET_LINK_LIBRARIES(test_quintic_interpolators_1d_nonuniform
                      sll_interpolators)

ADD_EXECUTABLE(test_odd_degree_interpolators_1d_nonuniform
               test_sll_odd_degree_spline_interpolator_1d_nonuniform.F90)

TARGET_LINK_LIBRARIES( test_odd_degree_interpolators_1d_nonuniform
                       sll_interpolators)

ADD_EXECUTABLE(test_odd_degree_interpolators_1d 
               test_sll_odd_degree_spline_interpolator_1d.F90)

TARGET_LINK_LIBRARIES( test_odd_degree_interpolators_1d
                       sll_interpolators)

ADD_EXECUTABLE(test_periodic_interpolators_1d test_periodic_interpolation.F90)

TARGET_LINK_LIBRARIES( test_periodic_interpolators_1d
                       sll_interpolators)
