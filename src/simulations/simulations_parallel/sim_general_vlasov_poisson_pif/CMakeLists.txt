ADD_LIBRARY( sll_general_vlasov_poisson_pif STATIC simulation_general_vlasov_poisson_pif.F90 )
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)

