ADD_LIBRARY(sll_qn_solvers STATIC sll_qns2d_with_finite_diff_seq.F90
                                  sll_qns2d_angular_spectral_method_seq.F90)

TARGET_LINK_LIBRARIES(sll_qn_solvers sll_fft sll_utilities sll_timer)

ADD_EXECUTABLE( test_qns2d_with_finite_diff_seq
                test_with_finite_diff_seq.F90)

TARGET_LINK_LIBRARIES( test_qns2d_with_finite_diff_seq sll_qn_solvers )

ADD_EXECUTABLE( test_qns2d_angular_spectral_method_seq 
                test_angular_spectral_method_seq.F90)

TARGET_LINK_LIBRARIES(test_qns2d_angular_spectral_method_seq sll_qn_solvers)

IF (MPI_MODULE_ENABLED)

   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   ADD_LIBRARY(sll_qn_solvers_parallel STATIC 
                  sll_qns2d_with_finite_diff_par.F90
                  sll_qns2d_angular_spectral_method_par.F90)

   TARGET_LINK_LIBRARIES(sll_qn_solvers_parallel sll_remap sll_fft)

   ADD_EXECUTABLE( test_qns2d_with_finite_diff_par 
                   test_with_finite_diff_par.F90)

   TARGET_LINK_LIBRARIES( test_qns2d_with_finite_diff_par sll_qn_solvers_parallel)

   ADD_EXECUTABLE( test_qns2d_angular_spectral_method_par 
                   test_angular_spectral_method_par.F90)

   TARGET_LINK_LIBRARIES( test_qns2d_angular_spectral_method_par 
                          sll_qn_solvers_parallel )

ENDIF (MPI_MODULE_ENABLED)

