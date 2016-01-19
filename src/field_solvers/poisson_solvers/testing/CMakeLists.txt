ADD_EXECUTABLE( test_poisson_2d_fft test_poisson_2d_fft.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_fft sll_poisson_solvers )
ADD_TEST( NAME poisson_2d_fft  COMMAND test_poisson_2d_fft     )
SET_TESTS_PROPERTIES( poisson_2d_fft PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_1d_fd test_poisson_1d_fd.F90 )
TARGET_LINK_LIBRARIES( test_poisson_1d_fd sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_1d_polar test_poisson_1d_polar.F90 )
TARGET_LINK_LIBRARIES( test_poisson_1d_polar sll_poisson_solvers )
ADD_TEST( NAME poisson_1d_polar COMMAND test_poisson_1d_polar_solver     )
SET_TESTS_PROPERTIES( poisson_1d_polar PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_2d_polar test_poisson_2d_polar.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_polar sll_poisson_solvers )
ADD_TEST( NAME poisson_2d_polar     COMMAND test_poisson_2d_polar     )
SET_TESTS_PROPERTIES( poisson_2d_polar PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_1d_periodic test_poisson_1d_periodic.F90 )
TARGET_LINK_LIBRARIES( test_poisson_1d_periodic sll_poisson_solvers )

# Cannot be compiled
#ADD_EXECUTABLE( test_poisson_1d_fem test_poisson_1d_fem.F90 )
#TARGET_LINK_LIBRARIES( test_poisson_1d_fem sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_2d_periodic test_poisson_2d_periodic.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_periodic sll_poisson_solvers )
ADD_TEST( NAME poisson_2d_periodic  COMMAND test_poisson_2d_periodic     )
SET_TESTS_PROPERTIES( poisson_2d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_2d_curvilinear test_poisson_2d_curvilinear.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_curvilinear sll_poisson_solvers )

ADD_EXECUTABLE( test_poisson_2d_polar_wrapper test_poisson_2d_polar_wrapper.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_polar_wrapper sll_poisson_solvers )
ADD_TEST( NAME poisson_2d_polar_wrapper     COMMAND test_poisson_2d_polar_wrapper     )
SET_TESTS_PROPERTIES( poisson_2d_polar_wrapper PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_3d_periodic test_poisson_3d_periodic.F90 )
TARGET_LINK_LIBRARIES( test_poisson_3d_periodic sll_poisson_solvers sll_fft )
ADD_TEST( NAME poisson_3d_periodic  COMMAND test_poisson_3d_periodic     )
SET_TESTS_PROPERTIES( poisson_3d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_tri_poisson test_tri_poisson.F90 )
TARGET_LINK_LIBRARIES( test_tri_poisson sll_poisson_solvers sll_mesh_calculus sll_file_io )
ADD_TEST( NAME tri_poisson  COMMAND test_tri_poisson     )
SET_TESTS_PROPERTIES( tri_poisson PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

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
ADD_TEST( NAME poisson_1d_periodic     COMMAND test_poisson_1d_periodic     )
SET_TESTS_PROPERTIES( poisson_1d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )



ADD_EXECUTABLE(test_fishpack test_fishpack.F90)
TARGET_LINK_LIBRARIES( test_fishpack sll_poisson_solvers)
ADD_TEST( NAME fishpack_interface COMMAND test_fishpack )
SET_TESTS_PROPERTIES( fishpack_interface PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

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
