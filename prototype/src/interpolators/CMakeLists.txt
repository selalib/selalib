IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

IF(STDF95)
  add_library(sll_interpolators STATIC
	sll_cubic_spline_interpolator_1d.F90
	sll_cubic_spline_interpolator_2d.F90
  )
ELSE()
  add_library(sll_interpolators STATIC
	sll_interpolators_1d_base.F90
	sll_cubic_spline_interpolator_1d.F90
	sll_finite_difference_interpolator_1d.F90
	sll_interpolators_2d_base.F90
	sll_cubic_spline_interpolator_2d.F90
  )
ENDIF()
add_dependencies(sll_interpolators
	sll_working_precision
	sll_memory
	sll_assertion
	sll_splines
)

add_executable(test_interpolators_2d unit_test.F90)
target_link_libraries(test_interpolators_2d
#	sll_weno
	sll_interpolators
	geometry_functions
	sll_splines
	util_constants
	sll_constants
	sll_assertion
	sll_memory
	sll_working_precision
)

#need to add test for finite difference 1D case.
add_executable(test_interpolators_1d unit_test_1d.F90)
target_link_libraries(test_interpolators_1d
#	sll_weno
	sll_interpolators
	geometry_functions
	sll_splines
	util_constants
	sll_constants
	sll_assertion
	sll_memory
	sll_working_precision
)

