#IF(FFT_LIB MATCHES "SLLFFT")

   ADD_LIBRARY(sll_qn_solvers STATIC sll_qn_solver_2d.F90)

   TARGET_LINK_LIBRARIES(sll_qn_solvers sll_fft sll_utilities sll_timer)

   ADD_EXECUTABLE( test_qn_solver_2d test_qn_solver_2d.F90)

   TARGET_LINK_LIBRARIES( test_qn_solver_2d sll_qn_solvers )

   ADD_TEST(NAME qns2d COMMAND test_qn_solver_2d)

   SET_TESTS_PROPERTIES(qns2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   IF (MPI_MODULE_ENABLED)

      SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})


      SET(QNS_PARALLEL_SRCS
        sll_qn_solver_2d.F90
        sll_qn_solver_2d_parallel.F90
        sll_qn_solver_3d_polar_parallel_x1.F90)

        SET(QNS_PARALLEL_SRCS ${QNS_PARALLEL_SRCS} 
          sll_qn_solver_3d_polar_parallel_x1_wrapper.F90)

      ADD_LIBRARY(sll_qn_solvers_parallel STATIC ${QNS_PARALLEL_SRCS})

      TARGET_LINK_LIBRARIES(sll_qn_solvers_parallel sll_remap sll_fft sll_poisson_solvers)
   
      ADD_EXECUTABLE( test_qn_solver_2d_parallel test_qn_solver_2d_parallel.F90)
   
      TARGET_LINK_LIBRARIES( test_qn_solver_2d_parallel sll_qn_solvers_parallel)

      ADD_EXECUTABLE( test_qn_solver_3d_polar_parallel_x1 
        unit_test_qn_solver_3d_polar_parallel_x1.F90)
   
      TARGET_LINK_LIBRARIES( test_qn_solver_3d_polar_parallel_x1 sll_qn_solvers_parallel)

      ADD_EXECUTABLE( test_qn_solver_3d_polar_parallel_x1_wrapper 
          unit_test_qn_solver_3d_polar_parallel_x1_wrapper.F90)
   
      TARGET_LINK_LIBRARIES( test_qn_solver_3d_polar_parallel_x1_wrapper 
          sll_qn_solvers_parallel)
        
   ENDIF (MPI_MODULE_ENABLED)

#ENDIF()

