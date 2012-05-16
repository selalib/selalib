add_library(geometry_functions STATIC
	../fields/geometry_functions.F90
)

add_library(interpolators STATIC
	sll_interpolators_1d_base.F90
	sll_cubic_spline_interpolator_1d.F90
	sll_interpolators_2d_base.F90
	sll_cubic_spline_interpolator_2d.F90
)

add_executable(test_interpolators_2d unit_test.F90)
target_link_libraries(test_interpolators_2d
	interpolators
	geometry_functions
	sll_splines
	util_constants
	sll_constants
	sll_assertion
	sll_memory)

add_executable(test_interpolators_1d unit_test_1d.F90)
target_link_libraries(test_interpolators_1d
#	sll_weno
	interpolators
	geometry_functions
	sll_splines
	util_constants
	sll_constants
	sll_assertion
	sll_memory)
