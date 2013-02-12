add_library (sll_lagrange_interpolation STATIC
		sll_lagrange_interpolation.F90)

add_dependencies(sll_lagrange_interpolation sll_memory 
                                            sll_assertion 
                                            sll_working_precision)

add_executable(test_lagrange_interpolation unit_test.F90)
target_link_libraries(test_lagrange_interpolation
			sll_lagrange_interpolation
                        sll_memory 
                        sll_assertion 
			sll_working_precision)
