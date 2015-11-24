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

ADD_EXECUTABLE( test_poisson_1d test_poisson_1d.F90 )
TARGET_LINK_LIBRARIES( test_poisson_1d sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_2d test_poisson_2d.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_2d_curvilinear test_poisson_2d_curvilinear.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_curvilinear sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_polar_fft test_poisson_polar_fft.F90 )
TARGET_LINK_LIBRARIES( test_poisson_polar_fft sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_3d_periodic_seq test_poisson_3d_periodic_seq.F90 )
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


ADD_EXECUTABLE(test_fishpack test_fishpack.F90)
TARGET_LINK_LIBRARIES( test_fishpack sll_poisson_solvers)

ADD_EXECUTABLE( test_poisson_2d_fem test_poisson_2d_fem.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_fem sll_poisson_solvers
  sll_memory
  sll_utilities
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES})

ADD_EXECUTABLE( test_hex_poisson test_hex_poisson.F90 )
TARGET_LINK_LIBRARIES( test_hex_poisson hex_poisson sll_utilities sll_constants )
