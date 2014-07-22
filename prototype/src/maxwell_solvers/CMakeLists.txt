IF(FFTW_FOUND AND FFTW_ENABLED)
  # this is very strange. Why compile fftw3.F90 here? this compilation should
  # depend on the configuration flags only!
   ADD_LIBRARY(sll_maxwell_solvers STATIC ../fft/fftw3.F90
                                          sll_maxwell_solvers_base.F90
                                          sll_maxwell_2d_fdtd.F90
                                          sll_maxwell_2d_pstd.F90
                                          sll_maxwell_2d_diga.F90
                                          sll_maxwell_3d_pstd.F90)

   TARGET_LINK_LIBRARIES(sll_maxwell_solvers sll_utilities ${FFTW_LIBRARIES})

   ADD_EXECUTABLE(test_maxwell_2d_pstd unit_test_2d_pstd.F90 )
   TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell_solvers)

   # ****************************************************************
   #                           PARALLEL MAXWELL SOLVERS
   # ***************************************************************

   ADD_LIBRARY(sll_maxwell_solvers STATIC sll_maxwell_solvers_base.F90 
                                          sll_maxwell_2d_fdtd.F90
                                          sll_maxwell_2d_diga.F90)

   TARGET_LINK_LIBRARIES(sll_maxwell_solvers sll_dg_fields 
                                             sll_boundary_condition_descriptors
                                             sll_utilities )

   ADD_EXECUTABLE(test_maxwell_2d_fdtd unit_test_2d_fdtd.F90)
   TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell_solvers)

   ADD_EXECUTABLE(test_maxwell_2d_diga_periodic unit_test_2d_diga_periodic.F90)
   TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_periodic sll_maxwell_solvers)

   ADD_EXECUTABLE(test_maxwell_2d_diga_wave unit_test_2d_diga_wave.F90)
   TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_wave sll_maxwell_solvers)

   ADD_EXECUTABLE(test_maxwell_2d_diga_discrete unit_test_2d_discrete.F90)
   TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_discrete sll_maxwell_solvers)

   IF(BUILD_TESTING)
      ADD_TEST(NAME maxwell_2d_dg_periodic COMMAND test_maxwell_2d_diga_periodic)
      SET_TESTS_PROPERTIES(maxwell_2d_dg_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
      ADD_TEST(NAME maxwell_2d_dg_silver_muller COMMAND test_maxwell_2d_diga_wave)
      SET_TESTS_PROPERTIES(maxwell_2d_dg_silver_muller PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDIF(BUILD_TESTING)


ENDIF(FFTW_FOUND AND FFTW_ENABLED)
