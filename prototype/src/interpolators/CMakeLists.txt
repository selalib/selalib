add_library(geometry_funcs STATIC
	../fields/geometry_functions.F90
)

add_library(interpolators STATIC
	../interpolators/sll_interpolators_2d_base.F90
	../interpolators/sll_cubic_spline_interpolator_2d.F90
)

#add_library(sll_interpolators_1d STATIC sll_interpolators_1d_base.F90)
#add_library(sll_interpolators_2d STATIC sll_interpolators_2d_base.F90)
#add_library(sll_cubic_spline_interpolator_1d STATIC sll_cubic_spline_interpolator_1d.F90)
#add_library(sll_cubic_spline_interpolator_2d STATIC sll_cubic_spline_interpolator_2d.F90)


#add_executable(test_interpolators_1d unit_test_1d.F90)
#target_link_libraries(test_interpolators_1d sll_interpolators_1d sll_WENO sll_splines sll_constants)
