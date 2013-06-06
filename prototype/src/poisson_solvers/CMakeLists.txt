SET(POISSON_LIBS sll_fft
                 sll_mapped_meshes
                 geometry_functions
                 sll_splines
                 sll_interpolators
                 sll_file_io
                 sll_utilities
                 sll_memory
                 sll_assertion
                 sll_timer
                 dfftpack)

SET(POISSON_SRCS sll_poisson_solvers.F90
                 sll_poisson_1d_periodic.F90
                 sll_poisson_2d_polar.F90)

IF(FFTW_FOUND AND FFTW_ENABLED)

   SET(POISSON_SRCS "${POISSON_SRCS};sll_poisson_2d_periodic_fftw.F90")

   ADD_EXECUTABLE(test_poisson_2d_fftw test_poisson_2d.F90
                                  sll_poisson_2d_periodic_fftw.F90)

   TARGET_LINK_LIBRARIES(test_poisson_2d_fftw ${POISSON_LIBS}
                                              ${FFTW_LIBRARIES}
                                              ${HDF5_LIBRARIES})
ELSE()

   SET(POISSON_SRCS "${POISSON_SRCS};sll_poisson_2d_periodic_fftpack.F90")

ENDIF()

IF(MUDPACK_ENABLED)

   SET(POISSON_SRCS ${POISSON_SRCS} sll_mudpack_base.F90
                                    sll_mudpack_cartesian.F90
                                    sll_mudpack_polar.F90
                                    sll_mudpack_colella.F90)

ENDIF(MUDPACK_ENABLED)
   
IF(FISHPACK_ENABLED)

   SET(POISSON_SRCS ${POISSON_SRCS} sll_fishpack.F90)

ENDIF(FISHPACK_ENABLED)

ADD_LIBRARY(sll_poisson_solvers STATIC ${POISSON_SRCS})

ADD_EXECUTABLE(test_poisson_1d test_poisson_1d.F90 sll_poisson_1d_periodic.F90)

ADD_LIBRARY(sll_poisson_2d_polar STATIC sll_poisson_2d_polar.F90)
ADD_EXECUTABLE(unit_test_polar unit_test_polar.F90 
               ../simulation/sll_polar_operators.F90)

ADD_EXECUTABLE(test_poisson_3d_periodic_seq test_poisson_3d_periodic_seq.F90
               sll_poisson_3d_periodic_seq.F90)

TARGET_LINK_LIBRARIES(unit_test_polar sll_poisson_2d_polar ${POISSON_LIBS})
TARGET_LINK_LIBRARIES(test_poisson_3d_periodic_seq ${POISSON_LIBS})
TARGET_LINK_LIBRARIES(test_poisson_1d ${POISSON_LIBS})

ADD_EXECUTABLE(test_poisson_2d_fem sll_poisson_2d_cart_fem.F90)

TARGET_LINK_LIBRARIES(test_poisson_2d_fem sll_memory 
                                          sll_assertion 
                                          ${LAPACK_LIBRARIES})

IF(FISHPACK_ENABLED)

   ADD_EXECUTABLE(test_fishpack test_fishpack.F90 sll_fishpack.F90)
   TARGET_LINK_LIBRARIES( test_fishpack sll_memory sll_assertion fishpack)

ENDIF()

IF(MUDPACK_ENABLED)

   MACRO(ADD_MUDPACK_PROGRAM TEST_NAME) 

      ADD_EXECUTABLE(test_mudpack_${TEST_NAME} test_mudpack_${TEST_NAME}.F90)
      
      TARGET_LINK_LIBRARIES(test_mudpack_${TEST_NAME} sll_memory 
                                                      sll_assertion 
                                                      sll_utilities
                                                      sll_file_io
                                                      sll_poisson_solvers
                                                      mudpack)

      ADD_TEST(NAME sll_mudpack_${TEST_NAME} COMMAND test_mudpack_${TEST_NAME} )
      SET(sll_mudpack_${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ENDMACRO(ADD_MUDPACK_PROGRAM)

   ADD_MUDPACK_PROGRAM(cartesian)
   ADD_MUDPACK_PROGRAM(polar)
   ADD_MUDPACK_PROGRAM(colella)


ENDIF()
# ****************************************************************************
#
#                           PARALLEL POISSON SOLVERS
#
# ****************************************************************************

IF(MPI_MODULE_ENABLED)

   ADD_LIBRARY(sll_poisson_parallel STATIC 
               sll_poisson_3d_periodic_par.F90 
               sll_poisson_2d_periodic_cartesian_par.F90 )

   ADD_DEPENDENCIES(sll_poisson_parallel sll_remap sll_collective)
  
   ADD_EXECUTABLE(test_poisson_3d_periodic_par test_poisson_3d_periodic_par.F90)

   TARGET_LINK_LIBRARIES( test_poisson_3d_periodic_par 
                          sll_poisson_parallel 
                          sll_remap 
                          sll_collective 
                          sll_fft 
                          sll_memory 
                          sll_timer 
                          sll_utilities 
                          sll_assertion)
ENDIF()

IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)

   # This test has hdf5 diagnostics. Maybe these should only be turned on
   # when debugging this specific module.
   ADD_EXECUTABLE( test_poisson_2d_per_cart_par 
                   test_poisson_2d_periodic_cartesian_par.F90 )

   TARGET_LINK_LIBRARIES(test_poisson_2d_per_cart_par
                         sll_poisson_parallel
                         sll_remap
                         sll_collective
                         sll_file_io_parallel
                         sll_fft
                         sll_memory
                         sll_timer
                         sll_utilities
                         sll_assertion )

ENDIF()
