add_library (sll_lagrange STATIC
		sll_lagrange_interpolation.F90)

add_executable(test_lagrange_interpolation unit_test_lagrange_interpolation.F90)
target_link_libraries(test_lagrange_interpolation
			sll_lagrange
			sll_working_precision)
