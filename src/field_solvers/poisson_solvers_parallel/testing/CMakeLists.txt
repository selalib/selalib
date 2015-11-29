SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

ADD_EXECUTABLE( test_poisson_3d_periodic_par test_poisson_3d_periodic_par.F90 )
TARGET_LINK_LIBRARIES( test_poisson_3d_periodic_par sll_poisson_parallel )

ADD_EXECUTABLE( test_poisson_2d_per_cart_par test_poisson_2d_periodic_cart_par.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_per_cart_par sll_poisson_parallel )

ADD_EXECUTABLE( test_poisson_polar_parallel test_poisson_polar_parallel.F90 )
TARGET_LINK_LIBRARIES( test_poisson_polar_parallel sll_poisson_parallel )

ADD_EXECUTABLE( test_poisson_2d_polar_parallel_solver test_poisson_2d_polar_parallel_solver.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_polar_parallel_solver sll_poisson_parallel )


IF(HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL)# AND FFT_LIB MATCHES "SLLFFT")
  SET( ARGS " " )
  ADD_MPI_TEST( poisson_3d_periodic_par test_poisson_3d_periodic_par ${PROCS} ${ARGS} ) 
  SET_TESTS_PROPERTIES( poisson_3d_periodic_par PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ADD_MPI_TEST( poisson_polar_parallel test_poisson_polar_parallel ${PROCS} ${ARGS} ) 
  SET_TESTS_PROPERTIES( poisson_polar_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ADD_MPI_TEST( poisson_2d_polar_parallel_solver test_poisson_2d_polar_parallel_solver ${PROCS} ${ARGS} ) 
  SET_TESTS_PROPERTIES( poisson_2d_polar_parallel_solver PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF()
IF(HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL AND FFT_LIB MATCHES "SLLFFT")
  ADD_MPI_TEST( poisson_per_cart_par_2d test_poisson_2d_per_cart_par ${PROCS} ${ARGS} ) 
  SET_TESTS_PROPERTIES( poisson_per_cart_par_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF()
