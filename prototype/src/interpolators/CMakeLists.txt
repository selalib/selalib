SET( INTERPOLATORS_SRC_FILES 
        sll_interpolators.F90
        sll_cubic_spline_interpolator_1d.F90
        sll_cubic_spline_interpolator_1d_nonuniform.F90
        sll_cubic_spline_interpolator_2d.F90
        sll_quintic_spline_interpolator_1d.F90
        sll_quintic_spline_interpolator_1d_nonuniform.F90
        sll_odd_degree_spline_interpolator_1d.F90
        sll_odd_degree_spline_interpolator_1d_nonuniform.F90
        sll_lagrange_interpolator_1d.F90 
        sll_periodic_interpolator_1d.F90 )

IF( NOT STDF95)

   SET( INTERPOLATORS_SRC_FILES ${INTERPOLATORS_SRC_FILES} 
        sll_interpolators_1d_base.F90
        sll_interpolators_2d_base.F90 )

ENDIF(NOT STDF95)


ADD_LIBRARY( sll_interpolators STATIC ${INTERPOLATORS_SRC_FILES} )

ADD_DEPENDENCIES( sll_interpolators
                  sll_splines
                  sll_periodic_interp
                  sll_fft
                  sll_working_precision
                  sll_memory
                  sll_assertion
                  sll_numerical_utilities 
                  sll_quintic_splines
                  sll_weno
                  sll_lagrange_interpolation )

IF(NOT STDF95)

   ADD_EXECUTABLE(test_interpolators_1d unit_test_1d.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../fields/geometry_functions.F90)

   TARGET_LINK_LIBRARIES( test_interpolators_1d
                          sll_interpolators
                          sll_splines
                          sll_fft
                          sll_utilities
                          sll_numerical_utilities
                          sll_constants
                          sll_assertion
                          sll_memory
                          sll_working_precision 
                          ${LAPACK_LIBRARIES} )
ENDIF()

ADD_EXECUTABLE(test_interpolators_2d unit_test.F90
               ${CMAKE_CURRENT_SOURCE_DIR}/../fields/geometry_functions.F90)

TARGET_LINK_LIBRARIES( test_interpolators_2d
                       sll_interpolators
                       sll_splines
                       sll_numerical_utilities  
                       sll_constants
                       sll_assertion
                       sll_memory
                       sll_working_precision )

ADD_EXECUTABLE(test_quintic_interpolators_1d 
               test_sll_quintic_spline_interpolator_1d.F90)

TARGET_LINK_LIBRARIES( test_quintic_interpolators_1d
                       sll_interpolators
                       sll_splines
                       sll_numerical_utilities
                       sll_constants
                       sll_assertion
                       sll_memory
                       sll_working_precision
                       ${LAPACK_LIBRARIES} )

ADD_EXECUTABLE(test_quintic_interpolators_1d_nonuniform 
               test_sll_quintic_spline_interpolator_1d_nonuniform.F90)

TARGET_LINK_LIBRARIES(test_quintic_interpolators_1d_nonuniform
                      sll_interpolators
                      sll_splines
                      sll_numerical_utilities
                      sll_constants
                      sll_assertion
                      sll_memory
                      sll_working_precision
                      ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(test_odd_degree_interpolators_1d_nonuniform
               test_sll_odd_degree_spline_interpolator_1d_nonuniform.F90)

TARGET_LINK_LIBRARIES( test_odd_degree_interpolators_1d_nonuniform
                       sll_interpolators
                       sll_splines
                       sll_numerical_utilities
                       sll_constants
                       sll_assertion
                       sll_memory
                       sll_working_precision
                       ${LAPACK_LIBRARIES} )

ADD_EXECUTABLE(test_odd_degree_interpolators_1d 
               test_sll_odd_degree_spline_interpolator_1d.F90)

TARGET_LINK_LIBRARIES( test_odd_degree_interpolators_1d
                       sll_interpolators
                       sll_splines
                       sll_numerical_utilities
                       sll_constants
                       sll_assertion
                       sll_memory
                       sll_working_precision
                       ${LAPACK_LIBRARIES} )

ADD_EXECUTABLE(test_periodic_interpolators_1d test_periodic_interpolation.F90)

TARGET_LINK_LIBRARIES( test_periodic_interpolators_1d
                       sll_interpolators
                       sll_periodic_interp
                       sll_splines
                       sll_numerical_utilities  
                       sll_fft
                       sll_utilities
                       sll_assertion
                       sll_memory
                       sll_working_precision
                       dfftpack
                       sll_lagrange_interpolation)
