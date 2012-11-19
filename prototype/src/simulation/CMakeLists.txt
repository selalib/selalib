#############################################################################################
###      Sequential simulations
#############################################################################################
# compile with openmp flag
IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic -qsmp=omp ${CMAKE_Fortran_FLAGS}")
ELSEIF(Fortran_COMPILER_NAME STREQUAL "gfortran")
  set(CMAKE_Fortran_FLAGS "-fopenmp ${CMAKE_Fortran_FLAGS}")
ELSEIF(Fortran_COMPILER_NAME STREQUAL "ifort")
  set(CMAKE_Fortran_FLAGS "-openmp ${CMAKE_Fortran_FLAGS}")
ENDIF()

#Here is defined the library containing the sequential simulations. 
#Please add your new simulation_MYMETHOD.F90 
############################################################

add_library(sll_simulation_1d STATIC simulation_1d.F90
            simulation_base_class.F90 
            simulation_VP1D_cartesian_non_unif.F90
            simulation_VP1D_curvilinear_analytic.F90
            simulation_vp1d_fourier_fem.F90)

#Add new dependencies concerning simulation_MYMETHOD.F90

add_dependencies( sll_simulation_1d
  cubic_non_uniform_splines
  sll_tridiagonal
  sll_splines
  sll_fft
  sll_constants
  sll_utilities
  sll_memory
  sll_assertion 	
  sll_working_precision
  )

#Here is defined the executables. 
#Please add your new test_MYMETHOD executable, built on unit_test_MYMETHOD.F90
#Add also the corresponding dependencies in a target_link_libraries() command
###############################################################

add_executable( test_VP1D_cartesian_non_unif unit_test_VP1D_cartesian_non_unif.F90 )

target_link_libraries( test_VP1D_cartesian_non_unif
  sll_simulation_1d
  cubic_non_uniform_splines
  sll_tridiagonal
  sll_constants
  sll_utilities
  sll_memory
  sll_assertion 	
  sll_working_precision
  )

add_executable( test_VP1D_curvilinear_analytic unit_test_VP1D_curvilinear_analytic.F90 )

target_link_libraries( test_VP1D_curvilinear_analytic
  sll_simulation_1d
  cubic_non_uniform_splines
  sll_tridiagonal
  sll_constants
  sll_utilities
  sll_memory
  sll_assertion 	
  sll_working_precision
  )




add_executable( test_vp1d_fourier_fem unit_test_vp1d_fourier_fem.F90 )

target_link_libraries( test_vp1d_fourier_fem
  sll_simulation_1d
  sll_splines
  sll_fft
  sll_constants
  sll_utilities
  sll_memory
  sll_assertion 	
  sll_working_precision
  )



#Commands to be updated in the above simulation class standard
#############################################################

#add_executable(test_simulation_1d user_geometry_functions.F90 unit_test.F90)
add_executable(polar_SL polar_SL.F90)
add_executable(curvilinear_SL curvilinear_SL.F90)

#target_link_libraries(test_simulation_1d
# 	sll_advection_field
#	sll_mapped_meshes
#	sll_interpolators
#	sll_splines
#	geometry_functions
#	sll_distribution_function
#	sll_scalar_field_2d	
#	sll_constants
#	sll_low_level_file_io
# 	sll_utilities
#	sll_memory
#	sll_assertion
#	sll_working_precision
#  )

target_link_libraries(polar_SL
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	geometry_functions
	sll_distribution_function
	sll_scalar_field_2d	
	sll_constants
	sll_low_level_file_io
  	sll_utilities
	sll_memory
	sll_assertion
  )
  
  target_link_libraries(curvilinear_SL
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	geometry_functions
	sll_distribution_function
	sll_scalar_field_2d	
	sll_constants
	sll_low_level_file_io
  	sll_utilities
	sll_memory
	sll_assertion
  )

add_executable(CG_polar
		CG_polar.F90
		#sll_poisson_2d_polar.F90
		sll_advection_polar.F90
		sll_polar_operators.F90
		)
add_dependencies(CG_polar sll_low_level_file_io)

FILE(COPY CG_data.txt DESTINATION ${CMAKE_BINARY_DIR})
target_link_libraries(CG_polar
	sll_splines
	sll_poisson_2d_polar
        sll_low_level_file_io
	sll_fft
	sll_tridiagonal
	sll_constants
	sll_memory
	sll_assertion
	sll_timer
	sll_utilities
	)
	

IF (0)
add_executable(Poisson_solver
		poisson_output.F90
		sll_poisson_2d_polar.F90
		)
target_link_libraries(Poisson_solver
	sll_fft
	sll_tridiagonal
	sll_constants
	sll_memory
	sll_assertion
	sll_utilities
	)
ENDIF()

IF(FFTPACK_ENABLED)
add_executable(VP1D_cart VP_uniform_cartesian_BSL.F90)
target_link_libraries(VP1D_cart
	sll_mapped_meshes
	sll_interpolators
	sll_fft
	sll_splines
	geometry_functions
	sll_distribution_function
	sll_scalar_field_2d	
	sll_constants
	sll_low_level_file_io
  	sll_utilities
	sll_memory
	sll_assertion 	
	sll_poisson_1d_periodic
	fftpack_module
	)

add_executable(VP1D_deltaf_cart VP1D_deltaf_BSL.F90)
target_link_libraries(VP1D_deltaf_cart
	sll_mapped_meshes
	sll_interpolators
	sll_fft
	sll_splines
	sll_periodic_interp
	geometry_functions
	sll_distribution_function
	sll_scalar_field_2d	
	sll_constants
	sll_low_level_file_io
  	sll_utilities
	sll_memory
	sll_assertion 	
	sll_poisson_1d_periodic
	fftpack_module
	${LAPACK_LIBRARIES}
	)

ENDIF()

#End of commands to be put in the class simulation standard
###################################################


#############################################################################################
###      Parallel simulations
#############################################################################################
# ATTENTION: From here below, the simulations/tests are parallel
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)

#Here is defined the library containing the parallel simulations. 
#Please add your new simulation_MYMETHOD.F90 
############################################################

add_library( sll_simulations STATIC simulation_base_class.F90 
  simulation_4d_vlasov_poisson_cartesian.F90 )

#Add new dependencies concerning simulation_MYMETHOD.F90

add_dependencies( sll_simulations
  sll_poisson_parallel
  sll_fft
  sll_electric_field_accumulators
  sll_scalar_field_4d
  sll_interpolators
  sll_splines
  sll_remap
  sll_collective
  sll_constants
  sll_low_level_file_io_parallel
  sll_utilities
  sll_memory
  sll_assertion 	
  sll_working_precision
  )

#Here is defined the executables. 
#Please add your new test_MYMETHOD executable, built on unit_test_MYMETHOD.F90
#Add also the corresponding dependencies in a target_link_libraries() command
###############################################################

add_executable( test_vp4d_sim unit_test_4d.F90 )

target_link_libraries( test_vp4d_sim
  sll_simulations
  sll_poisson_parallel
  sll_fft
  sll_electric_field_accumulators
  sll_scalar_field_4d
  sll_interpolators
  sll_splines
  sll_remap
  sll_collective
  sll_constants
  sll_low_level_file_io_parallel
  sll_utilities
  sll_memory
  sll_assertion 	
  sll_working_precision
  )

ENDIF() 
