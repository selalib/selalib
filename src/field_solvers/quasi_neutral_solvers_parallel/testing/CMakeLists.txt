# Parallel tests
IF(MPI_ENABLED AND HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL )
#IF(MPI_ENABLED AND HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL AND
#  FFT_LIB MATCHES "SLLFFT")
  
  ### BEGIN MACRO ###
  MACRO( QNS_PARALLEL_TEST TEST_NAME EXEC_NAME SOURCES )  
    SET( ARGS " " )
    SET( LIBS sll_quasi_neutral_solvers_parallel )
    ADD_EXECUTABLE( ${EXEC_NAME} ${SOURCES} )
    TARGET_LINK_LIBRARIES( ${EXEC_NAME} ${LIBS} )
    ADD_MPI_TEST( ${TEST_NAME} ${EXEC_NAME} ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES
      PASS_REGULAR_EXPRESSION "PASSED" )
  ENDMACRO( QNS_PARALLEL_TEST )
  ###  END MACRO  ###
  
  QNS_PARALLEL_TEST( qn_solver_2d_parallel
    test_qn_solver_2d_parallel
    test_qn_solver_2d_parallel.F90 )
 
  add_executable( test_qn_solver_2d_polar_par
    test_qn_solver_2d_polar_par.F90
    m_test_case_2d_dirichlet_1.F90
    m_test_case_2d_base.F90 )
  target_link_libraries( test_qn_solver_2d_polar_par
    sll_quasi_neutral_solvers_parallel )
  ADD_MPI_TEST( qn_solver_2d_polar_par
    test_qn_solver_2d_polar_par 1 " " )


# Test does not implement anything 
#  QNS_PARALLEL_TEST( qn_solver_3d_polar_parallel_x1
#    test_qn_solver_3d_polar_parallel_x1
#    test_qn_solver_3d_polar_parallel_x1.F90 )
  
  QNS_PARALLEL_TEST( qn_solver_3d_polar_parallel_x1_wrapper
    test_qn_solver_3d_polar_parallel_x1_wrapper
    test_qn_solver_3d_polar_parallel_x1_wrapper.F90 )
  
ENDIF()
