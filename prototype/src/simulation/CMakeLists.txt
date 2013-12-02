##############################################################################
###      Sequential simulations
##############################################################################
#Here is defined the library containing the sequential simulations. 
#Please add your new simulation_MYMETHOD.F90 
############################################################

ADD_LIBRARY(sll_simulation_1d STATIC simulation_1d.F90
                                     simulation_base_class.F90 
                                     simulation_VP1D_cartesian_non_unif.F90
                                     simulation_VP1D_curvilinear_analytic.F90
                                     simulation_vp1d_fourier_fem.F90
                                     simulation_2d_guiding_center_cartesian.F90
                                     simulation_2d_guiding_center_polar.F90)

TARGET_LINK_LIBRARIES( sll_simulation_1d sll_splines
                       sll_poisson_solvers
                       mudpack
                       sll_fft
		       sll_distribution_function
		       sll_general_coordinate_elliptic_solvers
		       deboor_splines
		       #sll_scalar_field_1d_alternative
		       sll_scalar_field_2d_alternative
                       sll_electric_field_accumulators
                       sll_scalar_field
                       sll_interpolators
                       sll_advection
                       sll_splines
		       sll_integration
		       sll_reduction
                       sll_constants
                       sll_coordinate_transformations
                       sll_logical_meshes
                       sll_utilities
               sll_common_array_initializers
		       sll_timer
                       sll_memory
                       sll_assertion 	
                       sll_working_precision)


# Define executables here. 
#Please add your new test_MYMETHOD executable, built on unit_test_MYMETHOD.F90
#Add also the corresponding dependencies in a target_link_libraries() command
###############################################################

MACRO(ADD_SIMULATION_1D TEST_NAME) 
   ADD_EXECUTABLE(test_${TEST_NAME} unit_test_${TEST_NAME}.F90)
   TARGET_LINK_LIBRARIES(test_${TEST_NAME} sll_simulation_1d)
ENDMACRO(ADD_SIMULATION_1D)

ADD_SIMULATION_1D(VP1D_cartesian_non_unif)
ADD_SIMULATION_1D(VP1D_curvilinear_analytic)
ADD_SIMULATION_1D(vp1d_fourier_fem)
ADD_SIMULATION_1D(2d_gc_cartesian)
ADD_SIMULATION_1D(2d_gc_polar)


#Commands to be updated in the above simulation class standard
#############################################################

ADD_EXECUTABLE(polar_SL polar_SL.F90)

TARGET_LINK_LIBRARIES(polar_SL sll_mapped_meshes
	                         sll_interpolators
	                         sll_constants
	                         sll_file_io
  	                         sll_utilities)

###

ADD_EXECUTABLE(curvilinear_SL curvilinear_SL.F90)
  
TARGET_LINK_LIBRARIES(curvilinear_SL sll_mapped_meshes
	                               sll_interpolators
	                               sll_constants
	                               sll_file_io
  	                               sll_utilities)
###

ADD_EXECUTABLE(CG_polar CG_polar.F90 sll_advection_polar.F90
		                         sll_polar_operators.F90)

TARGET_LINK_LIBRARIES(CG_polar sll_poisson_2d_polar
	                         sll_utilities
                               sll_file_io
	                         sll_fft
                               sll_splines
	                         sll_constants
	                         sll_timer)
###

ADD_EXECUTABLE(VP1D_cart VP_uniform_cartesian_BSL.F90)

TARGET_LINK_LIBRARIES(VP1D_cart sll_mapped_meshes
	                          sll_interpolators
	                          sll_fft
	                          sll_distribution_function
	                          sll_scalar_field	
	                          sll_constants
	                          sll_file_io
  	                          sll_utilities
	                          sll_poisson_solvers
	                          dfftpack)

###

#IF(CMAKE_BUILD_TYPE MATCHES "Release")
#OPENMP directives are activated only for release build

#OPENMP directives are activated only for release build
IF(CMAKE_BUILD_TYPE MATCHES "Release")

   ADD_EXECUTABLE(VP1D_deltaf_cart VP1D_deltaf_BSL.F90)
   TARGET_LINK_LIBRARIES(VP1D_deltaf_cart sll_mapped_meshes
	                                    sll_interpolators
	                                    sll_periodic_interp
	                                    sll_fft
	                                    sll_distribution_function
	                                    sll_scalar_field
                                          sll_constants
	                                    sll_file_io
  	                                    sll_utilities
	                                    sll_timer
	                                    sll_poisson_solvers
	                                    dfftpack)
   IF(HDF5_ENABLED)
     ADD_EXECUTABLE(VP1D_mm VP1D_mm.F90)
     TARGET_LINK_LIBRARIES(VP1D_mm sll_mapped_meshes
	                           sll_interpolators
	                           sll_periodic_interp
	                           sll_fft
	                           sll_distribution_function
	                           sll_scalar_field
	                           sll_constants
	                           sll_file_io
  	                           sll_utilities
	                           sll_timer
	                           sll_poisson_solvers
	                           dfftpack)
	ENDIF(HDF5_ENABLED)                           

ENDIF(CMAKE_BUILD_TYPE MATCHES "Release")

#End of commands to be put in the class simulation standard
###################################################

#############################################################################
###      Parallel simulations
#############################################################################
# ATTENTION: From here below, the simulations/tests are parallel
IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

   #Here is defined the library containing the parallel simulations. 
   #Please add there your new simulation_MYMETHOD.F90 
   ############################################################

   ADD_LIBRARY( sll_simulations STATIC 
     simulation_base_class.F90
     simulation_2d_vlasov_poisson_cartesian.F90
     simulation_4d_vlasov_poisson_cartesian.F90
     simulation_4d_vlasov_parallel_poisson_sequential_cartesian.F90
     simulation_4d_vlasov_poisson_general.F90
     simulation_4d_qns_general.F90
     #simulation_4d_mixed_qns_cartesian.F90
     simulation_4d_DK_hybrid.F90
     simulation_4d_vlasov_poisson_polar.F90
     simulation_4d_drift_kinetic_polar.F90
     simulation_4d_drift_kinetic_cartesian_finite_volume.F90
     simulation_6d_vlasov_poisson_cartesian.F90 )

   #Add new dependencies concerning simulation_MYMETHOD.F90
   SET(SIMULATION_DEPS sll_poisson_parallel
                       sll_poisson_solvers
                       sll_multigrid
                       sll_fft
		       sll_distribution_function
		       sll_general_coordinate_elliptic_solvers
		       deboor_splines
		       #sll_scalar_field_1d_alternative
		       sll_scalar_field_2d_alternative
                       sll_electric_field_accumulators
                       sll_scalar_field
                       sll_interpolators
                       sll_advection
                       sll_splines
                       sll_parallel_array_initializers	
                       sll_remap
                       sll_collective
                       sll_parallel_utilities
		       sll_integration
		       sll_reduction
                       sll_constants
                       sll_coordinate_transformations
                       sll_logical_meshes
                       sll_file_io_parallel
                       sll_utilities
		       sll_timer
                       sll_memory
                       sll_assertion 	
                       sll_working_precision)

   TARGET_LINK_LIBRARIES(sll_simulations ${SIMULATION_DEPS})

   MACRO(ADD_SIMULATION TEST_NAME) 
      ADD_EXECUTABLE(test_${TEST_NAME} unit_test_${TEST_NAME}.F90)
      TARGET_LINK_LIBRARIES(test_${TEST_NAME} sll_simulations)
   ENDMACRO(ADD_SIMULATION)

   #Define the executables/tests here.
   #Please add your new test_MYMETHOD executable, built on unit_test_MYMETHOD.F90
   #ADD_SIMULATION(MYMETHOD)
   ###############################################################

   ADD_SIMULATION(2d_vp_cartesian)
   ADD_SIMULATION(4d)
   ADD_SIMULATION(4d_vp_cartesian)
   ADD_SIMULATION(4d_vp_general)
   ADD_SIMULATION(4d_qns_general)
#   ADD_SIMULATION(4d_mixed_qns_cartesian)
   ADD_SIMULATION(4d_DK_hybrid)
   ADD_SIMULATION(4d_vp_polar)
   ADD_SIMULATION(4d_dk_polar)
   ADD_SIMULATION(4d_dk_cartesian)
   ADD_SIMULATION(6d)

   # Add postprocessing code for KEEN
   ADD_EXECUTABLE(VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f)
   TARGET_LINK_LIBRARIES(VPpostprocessing_drive_KEEN dfftpack)

ENDIF()
