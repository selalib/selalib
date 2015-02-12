ADD_LIBRARY(sll_time_splitting STATIC sll_time_splitting_base.F90
	                                   sll_const_coef_adv_2d.F90
	                                   sll_vlasov_poisson_adv_cart2d.F90
	                                   sll_time_splitting_coeff.F90)

TARGET_LINK_LIBRARIES(sll_time_splitting
                      sll_distribution_function
	                sll_interpolators
                      sll_utilities	
	                sll_file_io
                      sll_poisson_solvers)

IF(HDF5_ENABLED)

   ADD_EXECUTABLE(test_time_splitting unit_test.F90)
   TARGET_LINK_LIBRARIES(test_time_splitting sll_time_splitting)

ENDIF(HDF5_ENABLED)

ADD_EXECUTABLE(vp1d_bsl_time_split vp1d_bsl_time_split.F90)

TARGET_LINK_LIBRARIES(vp1d_bsl_time_split sll_time_splitting)
