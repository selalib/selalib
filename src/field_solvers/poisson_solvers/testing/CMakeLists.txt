# Save .mod files locally (this does not pollute the main "modules" directory)
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

ADD_EXECUTABLE( test_poisson_1d_polar test_poisson_1d_polar.F90 )
TARGET_LINK_LIBRARIES( test_poisson_1d_polar sll_poisson_solvers )
ADD_TEST( NAME poisson_1d_polar COMMAND test_poisson_1d_polar     )
SET_TESTS_PROPERTIES( poisson_1d_polar PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_2d_polar
  test_poisson_2d_polar.F90
  m_test_poisson_2d_polar_base.F90
  m_test_poisson_2d_polar_annulus_dirichlet.F90
  m_test_poisson_2d_polar_annulus_neumann_mode0.F90
  m_test_poisson_2d_polar_disk_dirichlet.F90 )

TARGET_LINK_LIBRARIES( test_poisson_2d_polar sll_poisson_solvers )
ADD_TEST( NAME poisson_2d_polar     COMMAND test_poisson_2d_polar     )
SET_TESTS_PROPERTIES( poisson_2d_polar PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_1d_periodic test_poisson_1d_periodic.F90 )
TARGET_LINK_LIBRARIES( test_poisson_1d_periodic sll_poisson_solvers )
ADD_TEST( NAME poisson_1d_periodic     COMMAND test_poisson_1d_periodic     )
SET_TESTS_PROPERTIES( poisson_1d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_1d_hmf test_poisson_1d_hmf.F90 )
TARGET_LINK_LIBRARIES( test_poisson_1d_hmf sll_poisson_solvers )
ADD_TEST( NAME poisson_1d_hmf     COMMAND test_poisson_1d_hmf     )
SET_TESTS_PROPERTIES( poisson_1d_hmf PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

# Cannot be compiled
#ADD_EXECUTABLE( test_poisson_1d_fem test_poisson_1d_fem.F90 )
#TARGET_LINK_LIBRARIES( test_poisson_1d_fem sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_2d_periodic test_poisson_2d_periodic.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_periodic sll_poisson_solvers )
ADD_TEST( NAME poisson_2d_periodic  COMMAND test_poisson_2d_periodic     )
SET_TESTS_PROPERTIES( poisson_2d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

IF(FFT_LIB MATCHES FFTW)
  ADD_EXECUTABLE( test_poisson_3d_periodic test_poisson_3d_periodic.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_3d_periodic sll_poisson_solvers sll_fft )
  ADD_TEST( NAME poisson_3d_periodic  COMMAND test_poisson_3d_periodic     )
  SET_TESTS_PROPERTIES( poisson_3d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(FFT_LIB MATCHES FFTW)

ADD_EXECUTABLE( test_poisson_2d_triangular test_poisson_2d_tri.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_triangular sll_poisson_solvers sll_mesh_calculus sll_file_io )
ADD_TEST( NAME poisson_2d_triangular  COMMAND test_poisson_2d_triangular     )
SET_TESTS_PROPERTIES( poisson_2d_triangular PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

IF(BUILD_MUDPACK)

  ADD_EXECUTABLE( test_poisson_2d_mudpack_cartesian test_mudpack_cartesian.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_mudpack_cartesian sll_poisson_solvers )
  ADD_TEST( NAME poisson_2d_mudpack_cartesian COMMAND test_poisson_2d_mudpack_cartesian )
  SET_TESTS_PROPERTIES( poisson_2d_mudpack_cartesian PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ADD_EXECUTABLE( test_poisson_2d_mudpack_polar test_mudpack_polar.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_mudpack_polar sll_poisson_solvers )
  ADD_TEST( NAME poisson_2d_mudpack_polar COMMAND test_poisson_2d_mudpack_polar )
  SET_TESTS_PROPERTIES( poisson_2d_mudpack_polar PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ADD_EXECUTABLE( test_poisson_2d_mudpack_curvilinear test_poisson_2d_mudpack_curvilinear.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_mudpack_curvilinear sll_poisson_solvers )
  ADD_TEST( NAME poisson_2d_mudpack_curvilinear COMMAND test_poisson_2d_mudpack_curvilinear )
  SET_TESTS_PROPERTIES( poisson_2d_mudpack_curvilinear PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(BUILD_MUDPACK)

ADD_EXECUTABLE( test_poisson_2d_fem test_poisson_2d_fem.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_fem sll_poisson_solvers
  sll_memory
  sll_utilities
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES})
ADD_TEST( NAME poisson_2d_fem COMMAND test_poisson_2d_fem )
SET_TESTS_PROPERTIES( poisson_2d_fem PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_hex_poisson test_hex_poisson.F90 )
TARGET_LINK_LIBRARIES( test_hex_poisson sll_poisson_solvers )
ADD_TEST( NAME hex_poisson COMMAND test_hex_poisson )
SET_TESTS_PROPERTIES( hex_poisson PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

IF(BUILD_MUDPACK)
  ADD_EXECUTABLE( test_poisson_2d_mudpack_curvilinear_old test_poisson_2d_mudpack_curvilinear_old.F90 )
  TARGET_LINK_LIBRARIES( test_poisson_2d_mudpack_curvilinear_old sll_poisson_solvers )
  ADD_TEST( NAME poisson_2d_mudpack_curvilinear_old COMMAND test_poisson_2d_mudpack_curvilinear_old )
  SET_TESTS_PROPERTIES( poisson_2d_mudpack_curvilinear_old PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_MUDPACK)

#-------------------------------------------------------------------------------
# Test FE Poisson solver based on smooth polar splines (dense version)
#-------------------------------------------------------------------------------
add_executable( test_poisson_2d_fem_sps_dense test_poisson_2d_fem_sps_dense.F90 )
target_link_libraries( test_poisson_2d_fem_sps_dense sll_poisson_solvers )

#-------------------------------------------------------------------------------
# Test FE Poisson solver based on smooth polar splines (stencil version)
#-------------------------------------------------------------------------------
add_executable( test_poisson_2d_fem_sps_stencil test_poisson_2d_fem_sps_stencil.F90 )
target_link_libraries( test_poisson_2d_fem_sps_stencil sll_poisson_solvers )
