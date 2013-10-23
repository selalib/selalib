ADD_LIBRARY(sll_distribution_function STATIC sll_distribution_function.F90
                                               sll_fdistribu4d_DK.F90
	                                       initial_functions.F90)

TARGET_LINK_LIBRARIES(sll_distribution_function sll_scalar_field
                                                sll_coordinate_transformations)

ADD_EXECUTABLE(test_distribution_function unit_test.F90)

TARGET_LINK_LIBRARIES(test_distribution_function sll_distribution_function)
