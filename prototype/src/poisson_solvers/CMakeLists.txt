
ADD_DEFINITIONS(-DMUDPACK)
SET(POISSON_SRCS sll_poisson_solvers.F90
                 sll_poisson_1d_base.F90
                 sll_poisson_1d_periodic_solver.F90
                 sll_poisson_1d_polar_solver.F90
                 sll_poisson_1d_periodic.F90
                 sll_poisson_2d_polar.F90
                 sll_poisson_2d_periodic_fftpack.F90
                 sll_poisson_2d_dirichlet_cartesian.F90
                 sll_poisson_2d_mudpack_solver.F90
                 #sll_poisson_2d_mudpack_curvilinear_solver_old.F90
                 sll_poisson_2d_base.F90
                 sll_poisson_2d_fft.F90
                 sll_poisson_2d_polar_solver.F90
                 sll_poisson_2d_periodic_solver.F90
                 sll_poisson_2d_elliptic_solver.F90
                 sll_poisson_3d_base.F90
                 sll_mudpack_base.F90
                 sll_mudpack_cartesian.F90
                 sll_mudpack_polar.F90
                 sll_mudpack_curvilinear.F90)

IF(FFTW_ENABLED AND FFTW_FOUND)
   SET(POISSON_SRCS 
       "${POISSON_SRCS};../fft/fftw3/fftw3.F90;sll_poisson_2d_periodic_fftw.F90")
ENDIF()

ADD_LIBRARY(sll_poisson_solvers STATIC ${POISSON_SRCS})

TARGET_LINK_LIBRARIES(sll_poisson_solvers sll_general_coordinate_elliptic_solvers 
                                          mudpack
                                          sll_memory 
                                          sll_assert 
                                          dfftpack )

IF(FFTW_ENABLED AND FFTW_FOUND)
   TARGET_LINK_LIBRARIES(sll_poisson_solvers ${FFTW_LIBRARIES})
ENDIF()


ADD_EXECUTABLE( test_poisson_2d_dirichlet_cartesian 
                test_poisson_2d_dirichlet_cartesian.F90 )

TARGET_LINK_LIBRARIES( test_poisson_2d_dirichlet_cartesian 
                       sll_poisson_solvers sll_utilities )

ADD_EXECUTABLE( test_poisson_2d_fft unit_test_poisson_2d_fft.F90)
TARGET_LINK_LIBRARIES( test_poisson_2d_fft sll_poisson_solvers)

ADD_EXECUTABLE(test_poisson_1d_polar_solver unit_test_poisson_1d_polar_solver.F90)
TARGET_LINK_LIBRARIES( test_poisson_1d_polar_solver sll_poisson_solvers)

ADD_EXECUTABLE( test_poisson_2d_polar_solver unit_test_poisson_2d_polar_solver.F90)
TARGET_LINK_LIBRARIES( test_poisson_2d_polar_solver sll_poisson_solvers)

ADD_EXECUTABLE(test_poisson_1d test_poisson_1d.F90 sll_poisson_1d_periodic.F90)
TARGET_LINK_LIBRARIES(test_poisson_1d sll_poisson_solvers)

ADD_EXECUTABLE(test_poisson_2d test_poisson_2d.F90)
TARGET_LINK_LIBRARIES(test_poisson_2d sll_poisson_solvers)

ADD_EXECUTABLE(test_poisson_polar_fft test_poisson_polar_fft.F90)
TARGET_LINK_LIBRARIES(test_poisson_polar_fft sll_poisson_solvers)

ADD_EXECUTABLE(test_poisson_3d_periodic_seq test_poisson_3d_periodic_seq.F90 
                                        sll_poisson_3d_periodic_seq.F90)
TARGET_LINK_LIBRARIES(test_poisson_3d_periodic_seq sll_fft)

ADD_EXECUTABLE(test_poisson_2d_fem sll_poisson_2d_fem.F90
                                   sll_poisson_2d_periodic_fem.F90
                                   test_poisson_2d_fem.F90)
TARGET_LINK_LIBRARIES(test_poisson_2d_fem sll_memory ${LAPACK_LIBRARIES})

IF(FISHPACK_ENABLED)

   ADD_EXECUTABLE(test_fishpack test_fishpack.F90 sll_fishpack.F90)
   TARGET_LINK_LIBRARIES( test_fishpack sll_memory sll_assert fishpack)

ENDIF()

### BEGIN MACRO

MACRO(ADD_MUDPACK_PROGRAM TEST_NAME) 

  ADD_EXECUTABLE(test_mudpack_${TEST_NAME} test_mudpack_${TEST_NAME}.F90)
    
  TARGET_LINK_LIBRARIES(test_mudpack_${TEST_NAME} ${POISSON_LIBS}
                                                sll_poisson_solvers 
                                                sll_file_io
                                                mudpack)

  ADD_TEST(NAME sll_mudpack_${TEST_NAME} COMMAND test_mudpack_${TEST_NAME} )
  SET(sll_mudpack_${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ENDMACRO(ADD_MUDPACK_PROGRAM)

### END MACRO

ADD_MUDPACK_PROGRAM(cartesian)
ADD_MUDPACK_PROGRAM(polar)
