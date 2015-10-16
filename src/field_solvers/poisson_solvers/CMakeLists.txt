ADD_DEFINITIONS(-DMUDPACK)
SET(POISSON_SRCS sll_m_poisson_1d_base.F90
                 sll_m_poisson_1d_periodic_solver.F90
                 sll_m_poisson_1d_polar_solver.F90
                 sll_m_poisson_1d_periodic.F90
                 sll_m_poisson_2d_polar.F90
                 sll_m_poisson_2d_periodic_fftpack.F90
                 sll_m_poisson_2d_dirichlet_cartesian.F90
                 sll_m_poisson_2d_mudpack.F90
                 sll_m_poisson_2d_mudpack_curvilinear_solver_old.F90
                 sll_m_poisson_2d_base.F90
                 sll_m_poisson_2d_fft.F90
                 sll_m_poisson_2d_polar_wrapper.F90
                 sll_m_poisson_2d_periodic_solver.F90
                 sll_m_poisson_2d_elliptic_solver.F90
                 sll_m_poisson_2d_curvilinear.F90
                 sll_m_poisson_3d_base.F90
                 sll_m_fishpack.F90
                 sll_m_mudpack.F90
                 sll_m_mudpack_curvilinear.F90
                 sll_m_poisson_1d_fourier.F90
                 sll_m_poisson_1d_fem.F90
                 sll_m_poisson_1d_fd.F90)

IF(FFTW_ENABLED AND FFTW_FOUND)
   SET(POISSON_SRCS ${POISSON_SRCS} sll_m_poisson_2d_periodic_fftw.F90)
ENDIF()

ADD_LIBRARY( sll_poisson_solvers STATIC ${POISSON_SRCS} )
TARGET_LINK_LIBRARIES( sll_poisson_solvers
  sll_general_coordinate_elliptic_solvers 
  sll_meshes
  sll_memory 
  sll_assert 
  mudpack
  fishpack
  dfftpack 
  ${LAPACK_LIBRARIES})

ADD_LIBRARY( hex_poisson STATIC sll_m_hex_poisson.F90 )
TARGET_LINK_LIBRARIES( hex_poisson sll_meshes sll_linear_solvers )

ADD_LIBRARY( tri_poisson STATIC sll_m_tri_poisson.F90 )
TARGET_LINK_LIBRARIES( tri_poisson sll_meshes )

# Tests
IF(BUILD_TESTING)

  IF(MPI_ENABLED)

  ADD_EXECUTABLE( test_poisson_2d_dirichlet_cartesian 
    test_poisson_2d_dirichlet_cartesian.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_dirichlet_cartesian 
    sll_poisson_solvers sll_utilities )

  ADD_EXECUTABLE( test_poisson_2d_fft test_poisson_2d_fft.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_fft sll_poisson_solvers )
  
  ADD_EXECUTABLE( test_poisson_1d_fd test_poisson_1d_fd.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_1d_fd sll_poisson_solvers )

  ADD_EXECUTABLE( test_poisson_1d_polar_solver test_poisson_1d_polar.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_1d_polar_solver sll_poisson_solvers )

  ADD_EXECUTABLE( test_poisson_2d_polar_solver test_poisson_2d_polar.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_polar_solver sll_poisson_solvers )

  ADD_EXECUTABLE( test_poisson_1d test_poisson_1d.F90 sll_m_poisson_1d_periodic.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_1d sll_poisson_solvers )

  ADD_EXECUTABLE( test_poisson_2d test_poisson_2d.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d sll_poisson_solvers )

  ADD_EXECUTABLE( test_poisson_2d_curvilinear test_poisson_2d_curvilinear.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_curvilinear sll_poisson_solvers )

  ADD_EXECUTABLE( test_poisson_polar_fft test_poisson_polar_fft.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_polar_fft sll_poisson_solvers )

  ADD_EXECUTABLE( test_poisson_3d_periodic_seq test_poisson_3d_periodic_seq.F90 
    sll_m_poisson_3d_periodic_seq.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_3d_periodic_seq sll_poisson_solvers sll_fft )




  ADD_EXECUTABLE( test_tri_poisson test_tri_poisson.F90 )
  TARGET_LINK_LIBRARIES( test_tri_poisson tri_poisson sll_mesh_calculus sll_file_io )

  ### BEGIN MACRO ###
  MACRO(ADD_MUDPACK_PROGRAM TEST_NAME)
    ADD_EXECUTABLE( test_mudpack_${TEST_NAME} test_mudpack_${TEST_NAME}.F90 )
    TARGET_LINK_LIBRARIES( test_mudpack_${TEST_NAME} mudpack
      sll_poisson_solvers sll_file_io )
    ADD_TEST( NAME sll_mudpack_${TEST_NAME} COMMAND test_mudpack_${TEST_NAME} )
    SET_TESTS_PROPERTIES( sll_mudpack_${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDMACRO(ADD_MUDPACK_PROGRAM)
  #### END MACRO ####

  # Ctests
  ADD_MUDPACK_PROGRAM( cartesian )
  ADD_MUDPACK_PROGRAM( polar )
  ADD_TEST( NAME poisson_solvers         COMMAND test_poisson_1d              )
  ADD_TEST( NAME poisson_3d_periodic_seq COMMAND test_poisson_3d_periodic_seq )
  SET_TESTS_PROPERTIES( poisson_3d_periodic_seq PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  IF(HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL AND FFT_LIB MATCHES "SLLFFT")
    SET( ARGS " " )
    ADD_MPI_TEST( poisson_3d_periodic_par test_poisson_3d_periodic_par ${PROCS} ${ARGS} ) 
    ADD_MPI_TEST( poisson_per_cart_par_2d test_poisson_2d_per_cart_par ${PROCS} ${ARGS} ) 
    SET_TESTS_PROPERTIES( poisson_3d_periodic_par PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
    SET_TESTS_PROPERTIES( poisson_per_cart_par_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF()

  ENDIF(MPI_ENABLED)

  ADD_EXECUTABLE(test_fishpack test_fishpack.F90)
  TARGET_LINK_LIBRARIES( test_fishpack sll_poisson_solvers) 
  
  ADD_EXECUTABLE( test_poisson_2d_fem sll_m_fem_2d.F90
    sll_m_fem_2d_periodic.F90 test_poisson_2d_fem.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_fem sll_memory 
                                           sll_utilities
                                           ${LAPACK_LIBRARIES} 
                                           ${BLAS_LIBRARIES})

  ADD_EXECUTABLE( test_hex_poisson test_hex_poisson.F90 )
  TARGET_LINK_LIBRARIES( test_hex_poisson hex_poisson sll_utilities sll_constants )

ENDIF(BUILD_TESTING)
