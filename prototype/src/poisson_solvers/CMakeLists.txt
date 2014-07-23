
SET(POISSON_SRCS sll_poisson_solvers.F90
                 sll_poisson_1d_periodic.F90
                 sll_poisson_2d_polar.F90
                 sll_poisson_2d_periodic_fftpack.F90)


IF(FFTW_ENABLED AND FFTW_FOUND)

   SET(POISSON_SRCS "${POISSON_SRCS};../fft/fftw3.F90;sll_poisson_2d_periodic_fftw.F90")

ENDIF()

IF(MUDPACK_ENABLED)
   ADD_DEFINITIONS(-DMUDPACK)
   SET(POISSON_SRCS ${POISSON_SRCS} sll_mudpack_base.F90
                                    sll_mudpack_cartesian.F90
                                    sll_mudpack_polar.F90
                                    sll_mudpack_curvilinear.F90)
ENDIF(MUDPACK_ENABLED)
   
IF(FISHPACK_ENABLED)

   SET(POISSON_SRCS ${POISSON_SRCS} sll_fishpack.F90)

ENDIF(FISHPACK_ENABLED)

SET(POISSON_SRCS ${POISSON_SRCS} sll_poisson_1d_base.F90
                                 sll_poisson_1d_periodic_solver.F90
                                 sll_poisson_1d_polar_solver.F90
                                 sll_poisson_2d_base.F90
                                 sll_poisson_2d_fft.F90
                                 sll_poisson_2d_polar_solver.F90
                                 sll_poisson_2d_periodic_solver.F90
                                 sll_poisson_2d_elliptic_solver.F90
                                 sll_poisson_3d_base.F90)

SET(POISSON_SRCS ${POISSON_SRCS} sll_poisson_2d_dirichlet_cartesian.F90 )
                                    
IF(MUDPACK_ENABLED)                                    
    SET(POISSON_SRCS ${POISSON_SRCS} sll_poisson_2d_mudpack_solver.F90
                                     sll_poisson_2d_mudpack_curvilinear_solver_old.F90)
ENDIF(MUDPACK_ENABLED)                                    
                                    
ADD_LIBRARY(sll_poisson_solvers STATIC ${POISSON_SRCS})

SET(POISSON_DEPS_LIBRARIES sll_general_coordinate_elliptic_solvers)

IF(MUDPACK_ENABLED)                                    
   SET(POISSON_DEPS_LIBRARIES ${POISSON_DEPS_LIBRARIES} mudpack) 
ENDIF(MUDPACK_ENABLED)                                    

IF(FFTW_ENABLED AND FFTW_FOUND)
   SET(POISSON_DEPS_LIBRARIES ${POISSON_DEPS_LIBRARIES} ${FFTW_LIBRARIES} )
ENDIF()

SET(POISSON_DEPS_LIBRARIES ${POISSON_DEPS_LIBRARIES} sll_memory sll_assertion dfftpack )

TARGET_LINK_LIBRARIES(sll_poisson_solvers ${POISSON_DEPS_LIBRARIES})

ADD_EXECUTABLE( test_poisson_2d_dirichlet_cartesian 
                test_poisson_2d_dirichlet_cartesian.F90 )

TARGET_LINK_LIBRARIES( test_poisson_2d_dirichlet_cartesian 
                       sll_poisson_solvers sll_utilities )

ADD_EXECUTABLE(
    test_poisson_2d_fft
    unit_test_poisson_2d_fft.F90)

TARGET_LINK_LIBRARIES(
    test_poisson_2d_fft
    sll_poisson_solvers)

ADD_EXECUTABLE(
    test_poisson_2d_polar_solver
    unit_test_poisson_2d_polar_solver.F90)

TARGET_LINK_LIBRARIES(
    test_poisson_2d_polar_solver
    sll_poisson_solvers
    sll_poisson_2d_polar)

ADD_EXECUTABLE(
    test_poisson_1d_polar_solver
    unit_test_poisson_1d_polar_solver.F90)

TARGET_LINK_LIBRARIES(
    test_poisson_1d_polar_solver
    sll_poisson_solvers)

SET(POISSON_LIBS sll_fft
                 sll_interpolators
                 sll_file_io
                 sll_utilities
                 sll_timer
		     sll_assertion
                 dfftpack)

ADD_EXECUTABLE(test_poisson_1d test_poisson_1d.F90 sll_poisson_1d_periodic.F90)
TARGET_LINK_LIBRARIES(test_poisson_1d ${POISSON_LIBS})

ADD_EXECUTABLE(test_poisson_2d test_poisson_2d.F90)

TARGET_LINK_LIBRARIES(test_poisson_2d sll_poisson_solvers)

ADD_LIBRARY(sll_poisson_2d_polar STATIC sll_poisson_2d_polar.F90)
TARGET_LINK_LIBRARIES(sll_poisson_2d_polar sll_fft sll_utilities)

ADD_EXECUTABLE(unit_test_polar unit_test_polar.F90 
               ../simulation/sll_polar_operators.F90)

TARGET_LINK_LIBRARIES(unit_test_polar sll_poisson_2d_polar ${POISSON_LIBS})

ADD_EXECUTABLE(test_poisson_polar_fft test_poisson_polar_fft.F90)
TARGET_LINK_LIBRARIES(test_poisson_polar_fft sll_poisson_2d_polar)

ADD_EXECUTABLE(test_poisson_3d_periodic_seq test_poisson_3d_periodic_seq.F90
               sll_poisson_3d_periodic_seq.F90)

TARGET_LINK_LIBRARIES(test_poisson_3d_periodic_seq ${POISSON_LIBS})

ADD_EXECUTABLE(test_poisson_2d_fem sll_poisson_2d_fem.F90
                                   sll_poisson_2d_periodic_fem.F90
                                   test_poisson_2d_fem.F90)

TARGET_LINK_LIBRARIES(test_poisson_2d_fem sll_memory ${LAPACK_LIBRARIES})

IF(FISHPACK_ENABLED)

   ADD_EXECUTABLE(test_fishpack test_fishpack.F90 sll_fishpack.F90)
   TARGET_LINK_LIBRARIES( test_fishpack sll_memory sll_assertion fishpack)

ENDIF()

IF(MUDPACK_ENABLED)

   MACRO(ADD_MUDPACK_PROGRAM TEST_NAME) 

      ADD_EXECUTABLE(test_mudpack_${TEST_NAME} test_mudpack_${TEST_NAME}.F90)
      
      TARGET_LINK_LIBRARIES(test_mudpack_${TEST_NAME} ${POISSON_LIBS}
                                                      sll_poisson_solvers 
                                                      sll_file_io
                                                      mudpack)

      ADD_TEST(NAME sll_mudpack_${TEST_NAME} COMMAND test_mudpack_${TEST_NAME} )
      SET(sll_mudpack_${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ENDMACRO(ADD_MUDPACK_PROGRAM)

   ADD_MUDPACK_PROGRAM(cartesian)
   ADD_MUDPACK_PROGRAM(polar)

ENDIF()
# ****************************************************************************
#
#                           PARALLEL POISSON SOLVERS
#
# ****************************************************************************

# This test has parallel diagnostics. Maybe these should only be turned on
# when debugging this specific module.
IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

  SET(POISSON_PARALLEL_SRCS
    sll_poisson_polar_parallel.F90 
    sll_poisson_3d_periodic_par.F90 
    sll_poisson_2d_periodic_cartesian_par.F90
    sll_quasi_neutral_3d_solver.F90)

    SET(POISSON_PARALLEL_SRCS ${POISSON_PARALLEL_SRCS} 
      sll_poisson_2d_polar_parallel_solver.F90)

   ADD_LIBRARY(sll_poisson_parallel STATIC ${POISSON_PARALLEL_SRCS})

   TARGET_LINK_LIBRARIES(sll_poisson_parallel sll_remap 
                                              sll_fft
                                              sll_file_io_parallel
                                              sll_collective
                                              sll_utilities
                                              sll_poisson_solvers)
  
   ADD_EXECUTABLE(test_poisson_3d_periodic_par test_poisson_3d_periodic_par.F90)
   TARGET_LINK_LIBRARIES(test_poisson_3d_periodic_par sll_poisson_parallel )

   ADD_EXECUTABLE(test_poisson_2d_per_cart_par 
                  test_poisson_2d_periodic_cartesian_par.F90 )
   TARGET_LINK_LIBRARIES(test_poisson_2d_per_cart_par sll_poisson_parallel)

   ADD_EXECUTABLE(test_poisson_polar_parallel 
                  test_poisson_polar_parallel.F90 )
   TARGET_LINK_LIBRARIES(test_poisson_polar_parallel sll_poisson_parallel)

   ADD_EXECUTABLE(test_poisson_polar_parallel_solver 
                  unit_test_poisson_polar_parallel_solver.F90 )
   TARGET_LINK_LIBRARIES(test_poisson_polar_parallel_solver sll_poisson_parallel)

   MESSAGE(STATUS "PASTIX_ENABLED:${PASTIX_ENABLED}")
   IF (PASTIX_ENABLED)
      ADD_EXECUTABLE(test_poisson_1d_pastix 
                     sll_poisson_1d_pastix.F90
                     test_poisson_1d_pastix.F90 )

      TARGET_LINK_LIBRARIES(test_poisson_1d_pastix sll_pastix)
   ENDIF (PASTIX_ENABLED)


ENDIF()
