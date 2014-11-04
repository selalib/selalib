##############################################################################
###      Sequential simulations
##############################################################################
#Here is defined the library containing the sequential simulations. 
#Please add your new simulation_MYMETHOD.F90 
############################################################

FILE(GLOB TXT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.txt")
FILE(COPY ${TXT_FILES} DESTINATION ${CMAKE_BINARY_DIR})
FILE(GLOB NML_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.nml")
FILE(COPY ${NML_FILES} DESTINATION ${CMAKE_BINARY_DIR})

ADD_DEFINITIONS(-DMUDPACK)

ADD_LIBRARY( sll_simulation_1d STATIC simulation_1d.F90 simulation_base_class.F90 )

TARGET_LINK_LIBRARIES( sll_simulation_1d 
                       sll_interpolators
                       sll_poisson_solvers
                       sll_reduction
                       sll_common_array_initializers
)

MACRO(ADD_SIMULATION_1D TEST_NAME SIM_NAME ) 
   ADD_EXECUTABLE(test_${TEST_NAME} simulation_${SIM_NAME}.F90
                                    unit_test_${TEST_NAME}.F90)
   TARGET_LINK_LIBRARIES(test_${TEST_NAME} sll_simulation_1d ${SIM_LIBS})

ENDMACRO(ADD_SIMULATION_1D)

#######################################################################
# Define executables here. 
# Please add your new test_MYMETHOD executable, built on 
# unit_test_MYMETHOD.F90 and simulation_MYMETHOD.F90
# Optionally add special dependencies by setting the DEP_LIBS variable
#######################################################################

ADD_SIMULATION_1D(VP1D_curvilinear_analytic VP1D_curvilinear_analytic)
ADD_SIMULATION_1D(vp1d_fourier_fem          vp1d_fourier_fem)
ADD_SIMULATION_1D(VP1D_cartesian_non_unif   VP1D_cartesian_non_unif)
IF(HDF5_ENABLED)
   SET(SIM_LIBS sll_advection)
ENDIF(HDF5_ENABLED)

ADD_EXECUTABLE(VP1D_deltaf_cart VP1D_deltaf_BSL.F90)
TARGET_LINK_LIBRARIES(VP1D_deltaf_cart sll_interpolators
	                                 sll_periodic_interpolation
	                                 sll_fft
	                                 sll_distribution_function
                                       sll_constants
	                                 sll_file_io
  	                                 sll_utilities
	                                 sll_timer
	                                 sll_poisson_solvers
	                                 dfftpack)
IF(HDF5_ENABLED)
   ADD_EXECUTABLE(VP1D_mm VP1D_mm.F90)
   TARGET_LINK_LIBRARIES(VP1D_mm sll_interpolators 
                                 sll_periodic_interpolation
                                 sll_fft 
                                 sll_distribution_function 
                                 sll_timer 
                                 sll_poisson_solvers 
                                 sll_constants 
                                 sll_file_io 
                                 sll_utilities 
                                 dfftpack)
ENDIF(HDF5_ENABLED)                           

#End of commands to be put in the class simulation standard
###################################################

#############################################################################
###      Parallel simulations
#############################################################################
# ATTENTION: From here below, the simulations/tests are parallel
IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

   #Here is defined the library containing the parallel simulations. 
   #Please add there your new simulation_MYMETHOD.F90 
   ############################################################

   ADD_LIBRARY( sll_simulations STATIC simulation_base_class.F90)

   #Add new dependencies concerning simulation_MYMETHOD.F90
   TARGET_LINK_LIBRARIES(sll_simulations 
                         sll_distribution_function_multipatch
                         #PN sll_general_coordinate_elliptic_solvers_multipatch
                         sll_scalar_field_multipatch
                         sll_parallel_array_initializers
                         sll_advection
                         sll_parallel_utilities
                         sll_poisson_solvers
                         sll_time_splitting
                         sll_file_io_parallel
                         sll_poisson_parallel
                         sll_fields
                         sll_fcisl
                         sll_derivative
                         sll_qn_solvers_parallel
                         sll_reduction
                         sll_gyroaverage
                         sll_pic_accumulators
                         sll_pic_utilities
                         sll_pic_initializers
                         sll_fft
                         sll_distribution_function
                         sll_general_coordinate_elliptic_solvers
                         sll_fields
                         sll_pic_particles
                         sll_random_deviate_generators
                         sll_field_accumulators
                         sll_interpolators
                         sll_point_to_point_communications	
                         sll_splines
                         sll_remap
                         sll_mesh_calculus
                         sll_integration
                         sll_constants
                         sll_coordinate_transformations
                         sll_meshes
                         sll_utilities
                         sll_timer
                         sll_memory
                         sll_assert 	
                         sll_working_precision
)


   MACRO(ADD_SIMULATION TEST_NAME SIM_NAME) 
      ADD_EXECUTABLE(test_${TEST_NAME} simulation_${SIM_NAME}.F90
                                       unit_test_${TEST_NAME}.F90)
      TARGET_LINK_LIBRARIES(test_${TEST_NAME} sll_simulations)
   ENDMACRO(ADD_SIMULATION)

   #Define the executables/tests here.
   #Please add your new test_MYMETHOD executable, built on unit_test_MYMETHOD.F90
   #and simulation_MYSIM.F90
   #ADD_SIMULATION(MYMETHOD MYSIM)

   ADD_SIMULATION(4d_qns_general    4d_qns_general)
   #PN ADD_SIMULATION(4d_qns_general_multipatch    4d_qns_general_multipatch)

   ADD_SIMULATION(4d_DK_hybrid              4d_DK_hybrid)
   ADD_SIMULATION(4d_dk_polar_multi_mu      4d_drift_kinetic_polar_multi_mu)
   ADD_SIMULATION(4d_dk_field_aligned_polar 4d_drift_kinetic_field_aligned_polar)

   ADD_SIMULATION(4d_vp_pic_cartesian       4d_vp_pic_cartesian)

#   ADD_SIMULATION(4d_mixed_qns_cartesian 4d_mixed_qns_cartesian)

   # Add postprocessing code for KEEN
   IF(Fortran_COMPILER MATCHES IBM)
      SET_SOURCE_FILES_PROPERTIES(VPpostprocessing_drive_KEEN.f 
                                  PROPERTIES COMPILE_FLAGS -qfixed)
   ENDIF(Fortran_COMPILER MATCHES IBM)

   ADD_EXECUTABLE(VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f)
   TARGET_LINK_LIBRARIES(VPpostprocessing_drive_KEEN dfftpack)

ENDIF()

ADD_EXECUTABLE(polar_SL polar_SL.F90)

TARGET_LINK_LIBRARIES(polar_SL sll_interpolators
	                       sll_constants
	                       sll_file_io
  	                       sll_utilities)

ADD_EXECUTABLE(curvilinear_SL curvilinear_SL.F90)
  
TARGET_LINK_LIBRARIES(curvilinear_SL sll_interpolators
	                             sll_constants
	                             sll_file_io
  	                             sll_utilities)

ADD_EXECUTABLE(CG_polar CG_polar.F90 sll_advection_polar.F90
		                         sll_polar_operators.F90)

TARGET_LINK_LIBRARIES(CG_polar sll_poisson_solvers
                               sll_utilities
                               sll_file_io
	                         sll_fft
                               sll_splines
	                         sll_constants
	                         sll_timer)

ADD_EXECUTABLE(VP1D_cart VP_uniform_cartesian_BSL.F90)

TARGET_LINK_LIBRARIES(VP1D_cart sll_interpolators
	                        sll_fft
	                        sll_distribution_function
	                        sll_constants
	                        sll_file_io
  	                        sll_utilities
	                        sll_poisson_solvers
	                        dfftpack)

ADD_EXECUTABLE(unit_test_polar unit_test_polar.F90 sll_polar_operators.F90)
TARGET_LINK_LIBRARIES(unit_test_polar sll_poisson_solvers)
