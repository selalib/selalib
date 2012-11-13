
IF(STDF95)
  add_library( sll_interpolators STATIC
	       sll_cubic_spline_interpolator_1d.F90
	       sll_cubic_spline_interpolator_2d.F90
               sll_quintic_spline_interpolator_1d.F90
             )
ELSE()
  add_library( sll_interpolators STATIC
	       sll_interpolators_1d_base.F90
	       sll_cubic_spline_interpolator_1d.F90
	       #sll_periodic_interpolator_1d.F90
	       #sll_finite_difference_interpolator_1d.F90
	       sll_interpolators_2d_base.F90
	       sll_cubic_spline_interpolator_2d.F90
               sll_quintic_spline_interpolator_1d.F90
               sll_odd_degree_spline_interpolator_1d.F90
             )
ENDIF()

add_dependencies( sll_interpolators
	          sll_working_precision
	          sll_memory
	          sll_assertion
	          sll_splines
                  sll_quintic_splines
                  sll_weno
		  #sll_periodic_interp
                )

add_executable(test_interpolators_1d unit_test_1d.F90
               ${CMAKE_SOURCE_DIR}/fields/geometry_functions.F90)

target_link_libraries( test_interpolators_1d
		       sll_interpolators
	               sll_splines
	               sll_constants
	               sll_assertion
	               sll_memory
	               sll_working_precision
                     )
add_executable(test_interpolators_2d unit_test.F90
               ${CMAKE_SOURCE_DIR}/fields/geometry_functions.F90)

target_link_libraries( test_interpolators_2d
	               sll_interpolators
	               sll_splines
	               sll_constants
	               sll_assertion
	               sll_memory
	               sll_working_precision
                     )

add_executable(test_quintic_interpolators_1d test_sll_quintic_spline_interpolator_1d.F90)
target_link_libraries( test_quintic_interpolators_1d
	               sll_interpolators
                       sll_splines
                       sll_numerical_utilities
	               sll_constants
	               sll_assertion
	               sll_memory
	               sll_working_precision
                     )

add_executable(test_odd_degree_interpolators_1d test_sll_odd_degree_spline_interpolator_1d.F90)
target_link_libraries( test_odd_degree_interpolators_1d
	               sll_interpolators
                       sll_splines
                       sll_numerical_utilities
	               sll_constants
	               sll_assertion
	               sll_memory
	               sll_working_precision
                       ${LAPACK_LIBRARIES}
                     )

#IF(FFTPACK_ENABLED)
#add_executable(test_periodic_interpolators_1d test_periodic_interpolation.F90)
#target_link_libraries( test_periodic_interpolators_1d
#                       sll_interpolators
#	               sll_periodic_interp
#                       sll_assertion
#	               sll_memory
#	               sll_working_precision
#		       fftpack_module
#                    )
#ENDIF()
