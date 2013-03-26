ADD_LIBRARY(sll_constants STATIC sll_constants.F90
	                              physical_constants.F90
	                              numeric_constants.F90)
ADD_DEPENDENCIES(sll_constants sll_working_precision)
ADD_EXECUTABLE(test_constants unit_test.F90)
TARGET_LINK_LIBRARIES(test_constants sll_constants sll_working_precision)
