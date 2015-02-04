IF(FFTW_FOUND AND FFTW_ENABLED)

  ADD_LIBRARY(sll_ampere_solvers STATIC sll_ampere_1d_pstd.F90)

  TARGET_LINK_LIBRARIES(sll_ampere_solvers 
                        sll_boundary_condition_descriptors
                        sll_utilities
                        sll_fft)

  IF(BUILD_TESTING)
     ADD_EXECUTABLE(test_ampere_1d_pstd unit_test_1d_pstd.F90 )
     TARGET_LINK_LIBRARIES(test_ampere_1d_pstd sll_ampere_solvers sll_file_io)
     ADD_TEST(NAME ampere_1d_pstd COMMAND test_ampere_1d_pstd)
     SET_TESTS_PROPERTIES(ampere_1d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  ENDIF()

ENDIF(FFTW_FOUND AND FFTW_ENABLED)
