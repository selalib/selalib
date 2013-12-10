ADD_LIBRARY(sll_linrood STATIC sll_linrood.F90)

TARGET_LINK_LIBRARIES(sll_linrood WENO_interp WENO_recon sll_scalar_field)

ADD_EXECUTABLE(test_linrood unit_test.F90)

TARGET_LINK_LIBRARIES(test_linrood sll_linrood sll_distribution_function)
