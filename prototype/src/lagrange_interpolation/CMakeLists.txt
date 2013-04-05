ADD_LIBRARY (sll_lagrange_interpolation STATIC
		sll_lagrange_interpolation.F90)

ADD_DEPENDENCIES(sll_lagrange_interpolation sll_memory 
                                            sll_assertion 
                                            sll_working_precision
					    sll_utilities)

ADD_EXECUTABLE(test_lagrange_interpolation unit_test.F90)
TARGET_LINK_LIBRARIES(test_lagrange_interpolation
			sll_lagrange_interpolation
                        sll_memory 
                        sll_assertion
			sll_working_precision
			sll_constants)
