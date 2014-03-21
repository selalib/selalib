SET( INTERPOLATORS_SRC_FILES 
        sll_cubic_spline_interpolator_1d.F90
        sll_cubic_spline_interpolator_1d_nonuniform.F90
        sll_quintic_spline_interpolator_1d.F90
        sll_quintic_spline_interpolator_1d_nonuniform.F90
        sll_odd_degree_spline_interpolator_1d.F90
        sll_odd_degree_spline_interpolator_1d_nonuniform.F90
        sll_lagrange_interpolator_1d.F90 
        sll_periodic_interpolator_1d.F90)

IF( NOT STDF95)

   SET( INTERPOLATORS_SRC_FILES ${INTERPOLATORS_SRC_FILES} 
        sll_interpolators_1d_base.F90
        sll_interpolators_2d_base.F90 
        sll_cubic_spline_interpolator_2d.F90
        sll_hermite_interpolator_2d.F90
	sll_arbitrary_degree_spline_interpolator_1d.F90
        sll_arbitrary_degree_spline_interpolator_2d.F90) 

ENDIF(NOT STDF95)

ADD_LIBRARY( sll_interpolators STATIC ${INTERPOLATORS_SRC_FILES} )

IF (NOT STDF95)
   TARGET_LINK_LIBRARIES( sll_interpolators sll_splines
                                            sll_periodic_interp
                                            sll_lagrange_interpolation
                                            sll_hermite_interpolation
                                            sll_fft
                                            sll_utilities 
   		                            deboor_splines 
                                            sll_coordinate_transformations
					    sll_timer
                                            ${LAPACK_LIBRARIES} )
ELSE()
   TARGET_LINK_LIBRARIES( sll_interpolators sll_splines
                                            sll_periodic_interp
                                            sll_lagrange_interpolation
                                            sll_hermite_interpolation
                                             sll_fft
                                            sll_utilities 
   		                            deboor_splines 
					    sll_file_io 
                                            ${LAPACK_LIBRARIES} )
ENDIF (NOT STDF95)


IF(NOT STDF95)

   ADD_EXECUTABLE(test_interpolators_1d unit_test_1d.F90)

   TARGET_LINK_LIBRARIES( test_interpolators_1d sll_interpolators)

   ADD_EXECUTABLE(test_arb_deg_spline_interpolators_2d 
                  unit_test_arbitrary_degree_2d.F90 )
   TARGET_LINK_LIBRARIES( test_arb_deg_spline_interpolators_2d
                          sll_interpolators
#                          sll_splines
#			    sll_file_io
#                          sll_utilities
#                          sll_constants
#                          sll_assertion
#                          sll_memory
#                          sll_working_precision 
#			  deboor_splines
)

   ADD_EXECUTABLE(test_arb_deg_spline_interpolators_1d 
                  unit_test_arbitrary_degree_1d.F90 )
   TARGET_LINK_LIBRARIES( test_arb_deg_spline_interpolators_1d
                          sll_interpolators sll_file_io)

   ADD_EXECUTABLE(test_interpolators_2d unit_test.F90)
   
   TARGET_LINK_LIBRARIES( test_interpolators_2d 
     sll_interpolators
     sll_utilities )

ENDIF()


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

IF(NOT STDF95)

ADD_EXECUTABLE(test_periodic_interpolators_1d test_periodic_interpolation.F90)

TARGET_LINK_LIBRARIES( test_periodic_interpolators_1d
                       sll_interpolators)
ENDIF(NOT STDF95)
