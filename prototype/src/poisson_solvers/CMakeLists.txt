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

ADD_EXECUTABLE(unit_test_polar unit_test_polar.F90 
                               sll_poisson_2d_polar.F90
	                            ../simulation/sll_polar_operators.F90)

TARGET_LINK_LIBRARIES(unit_test_polar ${POISSON_LIBS})

ADD_LIBRARY(sll_poisson_1d_periodic STATIC sll_poisson_1d_periodic.F90)
ADD_LIBRARY(sll_poisson_2d_polar    STATIC sll_poisson_2d_polar.F90)
ADD_LIBRARY(sll_poisson_3d_periodic STATIC sll_poisson_3d_periodic_seq.F90)

ADD_EXECUTABLE(test_poisson_3d_periodic_seq test_poisson_3d_periodic_seq.F90)

TARGET_LINK_LIBRARIES(test_poisson_3d_periodic_seq sll_poisson_3d_periodic
                                                   ${POISSON_LIBS})

ADD_EXECUTABLE(test_poisson_1d test_poisson_1d.F90 sll_poisson_1d_periodic.F90)

TARGET_LINK_LIBRARIES(test_poisson_1d ${POISSON_LIBS})

IF(FFTW_FOUND AND FFTW_ENABLED)


   ADD_LIBRARY(sll_poisson_solvers STATIC sll_poisson_1d_periodic.F90
                                          sll_poisson_2d_periodic_fftw.F90
                                          sll_poisson_2d_polar.F90)

   ADD_EXECUTABLE(test_poisson_2d test_poisson_2d.F90)

   TARGET_LINK_LIBRARIES(test_poisson_2d sll_poisson_solvers 
                                         ${POISSON_LIBS}
                                         ${FFTW_LIBRARIES})

ENDIF()

IF(MUDPACK_ENABLED)

   ADD_EXECUTABLE(test_mudpack test_mudpack.F90 sll_mudpack.F90)
   TARGET_LINK_LIBRARIES( test_mudpack sll_memory sll_assertion mudpack)

ENDIF()

IF(FISHPACK_ENABLED)

   ADD_EXECUTABLE(test_fishpack test_fishpack.F90 sll_fishpack.F90)
   TARGET_LINK_LIBRARIES( test_fishpack sll_memory sll_assertion fishpack)

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

ADD_EXECUTABLE(test_poisson_2d_fem sll_poisson_2d_cart_fem.F90)

TARGET_LINK_LIBRARIES(test_poisson_2d_fem sll_memory 
                                          sll_assertion 
                                          ${LAPACK_LIBRARIES})
