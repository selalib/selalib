add_library(sll_poisson_1d_periodic STATIC sll_poisson_1d_periodic.F90)

add_library(sll_poisson_2d_polar STATIC sll_poisson_2d_polar.F90)
target_link_libraries(sll_poisson_2d_polar
	sll_fft
	sll_tridiagonal
	sll_constants
	sll_memory
	sll_assertion
	sll_utilities
	)

add_executable(unit_test_polar
	unit_test_polar.F90
	../simulation/sll_polar_operators.F90
	#sll_poisson_2d_polar.F90
	)
target_link_libraries(unit_test_polar
	sll_fft
	sll_tridiagonal
	sll_constants
	sll_memory
	sll_assertion
	sll_utilities
	sll_poisson_2d_polar
	sll_splines
	)

MESSAGE(STATUS "FFTPACK_ENABLED:${FFTPACK_ENABLED}")

IF(FFTPACK_ENABLED)
  add_executable(test_poisson_1d test_poisson_1d.F90)
  target_link_libraries(test_poisson_1d
      sll_poisson_1d_periodic
      sll_scalar_field_1d
      sll_mapped_meshes
      geometry_functions
      sll_interpolators
      sll_low_level_file_io
      fftpack_module
      sll_utilities
      sll_assertion
      sll_memory
      sll_working_precision
  )

  add_definitions(-D_FFTPACK)
  add_library(sll_poisson_2d_periodic STATIC sll_poisson_2d_periodic.F90
	fft_module.F90)
  add_executable(test_poisson_2d test_poisson_2d.F90)
  target_link_libraries(test_poisson_2d 
                      sll_poisson_2d_periodic 
	              sll_scalar_field_2d
	              sll_constants
	              sll_working_precision
                      sll_assertion 
                      sll_memory 
                      sll_mapped_meshes
                      sll_scalar_field_1d
                      fftpack_module	 
                      geometry_functions	 
                      sll_low_level_file_io	
                      sll_utilities
  )
ELSE()

IF(FFTW_FOUND AND FFTW_ENABLED)
     ADD_LIBRARY(sll_poisson_2d_periodic STATIC sll_poisson_2d_periodic.F90)
     ADD_EXECUTABLE(test_poisson_2d test_poisson_2d.F90)
     TARGET_LINK_LIBRARIES(test_poisson_2d sll_poisson_2d_periodic 
                                           sll_mapped_meshes
                                   sll_constants
	                                   sll_working_precision
                                           sll_assertion 
                                           sll_memory 
                                           sll_low_level_file_io	
                                           sll_utilities
                                           sll_fft
                                           ${FFTW_LIBRARIES}
  )

ENDIF()

ENDIF()

IF(FISHPACK_ENABLED)
	add_executable(test_fishpack test_fishpack.F90 fishpack.F90)
	target_link_libraries( test_fishpack sll_memory sll_assertion fishpack)
ENDIF()

add_library(sll_poisson_3d_periodic_seq STATIC sll_poisson_3d_periodic_seq.F90)
add_executable(test_poisson_3d_periodic_seq test_poisson_3d_periodic_seq.F90)
target_link_libraries(test_poisson_3d_periodic_seq
	sll_poisson_3d_periodic_seq
	sll_fft
	sll_low_level_file_io	
        sll_memory
	sll_timer
	sll_utilities
	sll_assertion
)

# ****************************************************************************
#
#                           PARALLEL POISSON SOLVERS
#
# ****************************************************************************

#MPI_MODULE_ENABLED flag implies MPI_FOUND
IF(MPI_MODULE_ENABLED)

  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

  add_library(sll_poisson_parallel STATIC 
              sll_poisson_3d_periodic_par.F90 
              sll_poisson_2d_periodic_cartesian_par.F90 )

   ADD_DEPENDENCIES(sll_poisson_parallel sll_remap sll_collective)
  
IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)
  # This test has hdf5 diagnostics. Maybe these should only be turned on
  # when debugging this specific module.
  add_executable( test_poisson_2d_per_cart_par 
    test_poisson_2d_periodic_cartesian_par.F90 )

  target_link_libraries(test_poisson_2d_per_cart_par
    sll_poisson_parallel
    sll_remap
    sll_collective
    sll_low_level_file_io_parallel
    sll_fft
    sll_memory
    sll_timer
    sll_utilities
    sll_assertion )

ENDIF()

  add_executable(test_poisson_3d_periodic_par test_poisson_3d_periodic_par.F90)

  target_link_libraries( test_poisson_3d_periodic_par sll_poisson_parallel 
                         sll_remap sll_collective sll_fft sll_memory 
                         sll_timer sll_utilities sll_assertion)

ENDIF()


add_executable(test_poisson_2d_fem sll_poisson_2d_cart_fem.F90)

target_link_libraries(test_poisson_2d_fem sll_memory sll_assertion ${LAPACK_LIBRARIES})

