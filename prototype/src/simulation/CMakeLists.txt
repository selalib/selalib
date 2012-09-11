add_library(sll_simulation_1d STATIC simulation_1d.F90)

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
target_link_libraries(CG_polar
	sll_fft
	sll_tridiagonal
	sll_constants
	sll_splines
	sll_memory
	sll_assertion
	sll_timer
	sll_utilities
	sll_poisson_2d
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

ENDIF()
