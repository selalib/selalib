IF(FFTW_FOUND AND FFTW_ENABLED)

  ADD_LIBRARY(sll_ampere_solvers STATIC sll_ampere_vlasov_1d.F90)

  TARGET_LINK_LIBRARIES(sll_ampere_solvers 
                        sll_boundary_condition_descriptors
                        sll_utilities
                        sll_fft)

  IF(BUILD_TESTING)
     ADD_EXECUTABLE(test_ampere_vlasov_1d unit_test_ampere_vlasov_1d.F90 )
     TARGET_LINK_LIBRARIES(test_ampere_vlasov_1d sll_ampere_solvers 
                                                 sll_interpolators
                                                 sll_poisson_solvers
                                                 sll_file_io)
     ADD_TEST(NAME ampere_vlasov_1d COMMAND test_ampere_vlasov_1d)
     SET_TESTS_PROPERTIES(ampere_vlasov_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  ENDIF()

ENDIF(FFTW_FOUND AND FFTW_ENABLED)
