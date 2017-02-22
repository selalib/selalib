SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

IF (PGI)
  ADD_LIBRARY(sll_maxwell_parallel STATIC 
              sll_m_maxwell_2d_periodic_cartesian_par_pgi.F90 )
ELSE(PGI)
  ADD_LIBRARY(sll_maxwell_parallel STATIC 
              sll_m_maxwell_2d_periodic_cartesian_par.F90 )
ENDIF(PGI)

TARGET_LINK_LIBRARIES(sll_maxwell_parallel sll_remap sll_fft )
  
IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()
