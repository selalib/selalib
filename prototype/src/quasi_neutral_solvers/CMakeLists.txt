IF(FFT_LIB MATCHES "SLLFFT")

   ADD_LIBRARY(sll_qn_solvers STATIC sll_qn_solver_2d.F90)

   TARGET_LINK_LIBRARIES(sll_qn_solvers sll_fft sll_utilities sll_timer)

   ADD_EXECUTABLE( test_qn_solver_2d test_qn_solver_2d.F90)

   TARGET_LINK_LIBRARIES( test_qn_solver_2d sll_qn_solvers )

   ADD_TEST(NAME qns2d COMMAND test_qn_solver_2d)

   SET_TESTS_PROPERTIES(qns2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   IF (MPI_MODULE_ENABLED)

      SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

      ADD_LIBRARY(sll_qn_solvers_parallel STATIC 
                     sll_qn_solver_2d.F90 sll_qn_solver_2d_parallel.F90)

      TARGET_LINK_LIBRARIES(sll_qn_solvers_parallel sll_remap sll_fft)
   
      ADD_EXECUTABLE( test_qn_solver_2d_parallel test_qn_solver_2d_parallel.F90)
   
      TARGET_LINK_LIBRARIES( test_qn_solver_2d_parallel sll_qn_solvers_parallel)

   ENDIF (MPI_MODULE_ENABLED)

ENDIF()

