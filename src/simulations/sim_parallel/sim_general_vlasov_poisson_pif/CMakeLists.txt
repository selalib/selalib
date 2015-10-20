
ADD_LIBRARY( sll_general_vlasov_poisson_pif STATIC simulation_general_vlasov_poisson_pif.F90 )
TARGET_LINK_LIBRARIES (sll_general_vlasov_poisson_pif sll_sim_parallel
					      sll_pif
					      sll_timer
					      sll_memory
					      sll_utilities
					      sll_constants
					      sll_sobol
					      sll_prob
					      sll_pic_visu
					      sll_pic_visu_parallel
					      sll_collective
					      sll_moment_matching
					      sll_pic_utilities
					      sll_descriptors
					      sll_splitting_methods)
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
