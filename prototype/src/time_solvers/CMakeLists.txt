ADD_LIBRARY(sll_time_splitting STATIC sll_time_splitting_base.F90
	                                   sll_const_coef_adv_2d.F90
	                                   sll_vlasov_poisson_adv_cart2d.F90)

ADD_DEPENDENCIES(sll_time_splitting
                 sll_distribution_function
	              sll_interpolators
	              sll_splines
                 sll_numerical_utilities	
	              sll_low_level_file_io
	              sll_memory
	              sll_assertion
	              sll_working_precision
                 sll_poisson_1d_periodic)

ADD_EXECUTABLE(test_time_splitting unit_test.F90)

TARGET_LINK_LIBRARIES(test_time_splitting sll_time_splitting
	                                       sll_interpolators
	                                       sll_low_level_file_io
	                                       sll_splines
                                          sll_numerical_utilities	
	                                       sll_memory
	                                       sll_assertion
	                                       sll_working_precision)

IF(FFTPACK_ENABLED)

  ADD_EXECUTABLE(vp1d_bsl_time_split vp1d_bsl_time_split.F90)

  TARGET_LINK_LIBRARIES(vp1d_bsl_time_split sll_mapped_meshes
                                            sll_interpolators
                                            sll_splines
                                            geometry_functions
                                            sll_distribution_function
                                            sll_scalar_field_2d
                                            sll_time_splitting
                                            sll_numerical_utilities	
                                            sll_constants
                                            sll_low_level_file_io
                                            sll_memory
                                            sll_assertion
                                            sll_working_precision
                                            sll_utilities
                                            sll_poisson_1d_periodic
                                            dfftpack)


ENDIF()
