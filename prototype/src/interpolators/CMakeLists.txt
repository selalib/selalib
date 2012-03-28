add_library(sll_interpolators STATIC
	sll_interpolators_base.F90
	sll_cubic_spline_interpolator_2d.F90
)

add_library(geometry_funcs STATIC
	../fields/geometry_functions.F90
)

add_executable(test_interpolators unit_test.F90)
target_link_libraries(test_interpolators sll_interpolators sll_mesh_types sll_splines sll_memory sll_constants sll_assertion)
