##############################################################################
###      Sequential simulations
##############################################################################
#Here is defined the library containing the sequential simulations. 
#Please add your new simulation_MYMETHOD.F90 
############################################################

IF (STDF95)
add_library(sll_simulation_1d STATIC simulation_1d.F90
            simulation_VP1D_cartesian_non_unif.F90
            simulation_VP1D_curvilinear_analytic.F90
            simulation_vp1d_fourier_fem.F90)
ELSE()
add_library(sll_simulation_1d STATIC simulation_1d.F90
            simulation_base_class.F90 
            simulation_VP1D_cartesian_non_unif.F90
            simulation_VP1D_curvilinear_analytic.F90
            simulation_vp1d_fourier_fem.F90)
ENDIF (STDF95)

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

# Define executables here. 
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
#	sll_file_io
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
#	sll_distribution_function
#	sll_scalar_field_2d	
	sll_constants
	sll_file_io
  	sll_utilities
	sll_memory
	sll_assertion
  )
  
  target_link_libraries(curvilinear_SL
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	geometry_functions
#	sll_distribution_function
#	sll_scalar_field_2d	
	sll_constants
	sll_file_io
  	sll_utilities
	sll_memory
	sll_assertion
  )

add_executable(CG_polar CG_polar.F90
		                  sll_advection_polar.F90
		                  sll_polar_operators.F90)

add_dependencies(CG_polar sll_file_io sll_fft sll_poisson_2d_polar)

#FILE(COPY CG_data.txt DESTINATION ${CMAKE_BINARY_DIR})

target_link_libraries(CG_polar 
	                            sll_poisson_2d_polar
	                            sll_utilities
                               sll_file_io
	                            sll_fft
                               sll_splines
	                            sll_constants
	                            sll_memory
	                            sll_assertion
	                            sll_timer
                                    ${HDF5_LIBRARIES}
)

IF (0)
add_executable(Poisson_solver
		poisson_output.F90
		sll_poisson_2d_polar.F90
		)
target_link_libraries(Poisson_solver
	sll_fft
#	sll_tridiagonal
	sll_constants
	sll_memory
	sll_assertion
	sll_utilities
	)
ENDIF()

IF(NOT STDF95)
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
	sll_file_io
  	sll_utilities
	sll_memory
	sll_assertion 	
	sll_poisson_solvers
	dfftpack
                                    ${HDF5_LIBRARIES}
	)


IF(CMAKE_BUILD_TYPE MATCHES "Release")

add_executable(VP1D_deltaf_cart VP1D_deltaf_BSL.F90)
target_link_libraries(VP1D_deltaf_cart
	sll_mapped_meshes
	sll_interpolators
	sll_periodic_interp
	sll_splines
	sll_fft
	geometry_functions
	sll_distribution_function
	sll_scalar_field_2d	
	sll_constants
	sll_file_io
  	sll_utilities
	sll_timer
	sll_memory
	sll_assertion 	
	sll_poisson_solvers
	dfftpack
	${LAPACK_LIBRARIES}
	)

add_executable(VP1D_mm VP1D_mm.F90)
target_link_libraries(VP1D_mm
	sll_mapped_meshes
	sll_interpolators
	sll_periodic_interp
	sll_splines
	sll_fft
	geometry_functions
	sll_distribution_function
	sll_scalar_field_2d	
	sll_constants
	sll_file_io
  	sll_utilities
	sll_timer
	sll_memory
	sll_assertion 	
	sll_poisson_solvers
	dfftpack
	${LAPACK_LIBRARIES}
	)

ENDIF(CMAKE_BUILD_TYPE MATCHES "Release")

ENDIF(NOT STDF95)

#End of commands to be put in the class simulation standard
###################################################


#############################################################################
###      Parallel simulations
#############################################################################
# ATTENTION: From here below, the simulations/tests are parallel
IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)

   #Here is defined the library containing the parallel simulations. 
   #Please add there your new simulation_MYMETHOD.F90 
   ############################################################

   ADD_LIBRARY( sll_simulations STATIC 
     simulation_base_class.F90
     simulation_4d_vlasov_poisson_cartesian.F90
     simulation_4d_vlasov_poisson_general.F90
     simulation_4d_qns_general.F90
     simulation_4d_vlasov_poisson_polar.F90
     simulation_4d_drift_kinetic_cartesian_finite_volume.F90
     simulation_6d_vlasov_poisson_cartesian.F90)

   #Add new dependencies concerning simulation_MYMETHOD.F90
   SET(SIMULATION_DEPS sll_poisson_parallel
                       sll_fft
		       sll_general_coordinate_quasi_neutral_solvers
		       sll_scalar_field_2d_alternative
                       sll_electric_field_accumulators
                       sll_scalar_field
                       sll_interpolators
                       sll_splines
                       sll_utilities
                       sll_parallel_array_initializers	
                       sll_remap
                       sll_collective
		       sll_integration
                       sll_constants
                       sll_utilities
                       sll_coordinate_transformations
                       sll_logical_meshes
                       sll_file_io_parallel
                       sll_memory
                       sll_assertion 	
                       sll_working_precision
		       deboor_splines)

   ADD_DEPENDENCIES(sll_simulations ${SIMULATION_DEPS})

   MACRO(ADD_SIMULATION TEST_NAME) 
      ADD_EXECUTABLE(test_${TEST_NAME} unit_test_${TEST_NAME}.F90)
      TARGET_LINK_LIBRARIES(test_${TEST_NAME} sll_simulations ${SIMULATION_DEPS})
   ENDMACRO(ADD_SIMULATION)

   #Define the executables/tests here.
   #Please add your new test_MYMETHOD executable, built on unit_test_MYMETHOD.F90
   #ADD_SIMULATION(MYMETHOD)
   ###############################################################

   ADD_SIMULATION(4d)
   ADD_SIMULATION(4d_vp_general)
   ADD_SIMULATION(4d_qns_general)
   ADD_SIMULATION(4d_vp_polar)
   ADD_SIMULATION(4d_dk_cartesian)
   ADD_SIMULATION(6d)

ENDIF()
