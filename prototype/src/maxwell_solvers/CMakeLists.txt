IF(FFTW_FOUND AND NOT STDF95)

   ADD_LIBRARY(sll_maxwell STATIC sll_maxwell.F90
                                  sll_maxwell_2d_fdtd.F90
                                  sll_maxwell_2d_pstd.F90)

   ADD_EXECUTABLE(test_maxwell_2d_pstd unit_test_2d_pstd.F90 )
   
   TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell
                                              sll_assertion 
                                              sll_constants 
                                              sll_memory 
                                              sll_collective 
                                              sll_utilities 
                                              ${FFTW_LIBRARIES})
ELSE()

   ADD_LIBRARY(sll_maxwell STATIC sll_maxwell.F90 sll_maxwell_2d_fdtd.F90)

ENDIF(FFTW_FOUND AND NOT STDF95)

ADD_EXECUTABLE(test_maxwell_2d_fdtd unit_test_2d_fdtd.F90)

ADD_DEPENDENCIES(sll_maxwell sll_assertion 
                             sll_constants 
                             sll_memory 
                             sll_collective 
                             sll_utilities )

TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell
                                           sll_assertion 
                                           sll_constants 
                                           sll_memory 
                                           sll_collective 
                                           sll_utilities )


# ****************************************************************************
#
#                           PARALLEL MAXWELL SOLVERS
#
# ****************************************************************************

IF(MPI_MODULE_ENABLED AND HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

  ADD_LIBRARY(sll_maxwell_parallel STATIC 
              sll_maxwell_2d_periodic_cartesian_par.F90 )

  ADD_DEPENDENCIES(sll_maxwell_parallel sll_remap sll_collective)
  
  ADD_EXECUTABLE(test_maxwell_2d_per_cart_par test_maxwell_2d_periodic_cartesian_par.F90 )

  TARGET_LINK_LIBRARIES(test_maxwell_2d_per_cart_par sll_maxwell_parallel
                                                     sll_remap
                                                     sll_collective
                                                     sll_low_level_file_io_parallel
                                                     sll_fft
                                                     sll_memory
                                                     sll_timer
                                                     sll_utilities
                                                     sll_assertion )
ENDIF()
