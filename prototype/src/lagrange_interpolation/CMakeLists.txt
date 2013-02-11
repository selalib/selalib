add_library (sll_lagrange_interpolation STATIC
		sll_lagrange_interpolation.F90)

add_executable(test_lagrange_interpolation unit_test.F90)
target_link_libraries(test_lagrange_interpolation
			sll_lagrange_interpolation
			sll_working_precision)


#test with function in lagrange_interpolation_test
# add_executable(test_lagrange_interpolation unit_test_lagrange_interpolation.F90)
# target_link_libraries(test_lagrange_interpolation
# 			sll_lagrange_interpolation
# 			sll_working_precision)
# 
# add_executable(test_lagrange_interpolation_2d unit_test_lagrange_interpolation_2d.F90)
# target_link_libraries(test_lagrange_interpolation_2d
# 			sll_lagrange_interpolation
# 			sll_working_precision)
# 
# add_executable(test_lagrange_interpolation_3d unit_test_lagrange_interpolation_3d.F90)
# target_link_libraries(test_lagrange_interpolation_3d
# 			sll_lagrange_interpolation
# 			sll_working_precision)
