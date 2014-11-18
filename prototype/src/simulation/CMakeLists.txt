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

IF(HDF5_ENABLED)
   SET(SIM_LIBS sll_advection)
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

#YG   ADD_SIMULATION(4d_qns_general    4d_qns_general)
#PN   ADD_SIMULATION(4d_qns_general_multipatch    4d_qns_general_multipatch)

#PN   ADD_SIMULATION(4d_DK_hybrid              4d_DK_hybrid)
#PN   ADD_SIMULATION(4d_dk_polar               4d_drift_kinetic_polar)
#PN   ADD_SIMULATION(4d_dk_polar_one_mu        4d_drift_kinetic_polar_one_mu)
#PN   ADD_SIMULATION(4d_dk_polar_multi_mu      4d_drift_kinetic_polar_multi_mu)
#PN   ADD_SIMULATION(4d_dk_field_aligned_polar 4d_drift_kinetic_field_aligned_polar)
#PN   ADD_SIMULATION(6d                        6d_vlasov_poisson_cartesian)

#PN   ADD_SIMULATION(4d_vp_pic_cartesian       4d_vp_pic_cartesian)

#SH   ADD_SIMULATION(2d_guiding_center_pic_cartesian 2d_guiding_center_pic_cartesian)


#   ADD_SIMULATION(4d_mixed_qns_cartesian 4d_mixed_qns_cartesian)

#YG   # Add postprocessing code for KEEN
#YG   IF(Fortran_COMPILER MATCHES IBM)
#YG      SET_SOURCE_FILES_PROPERTIES(VPpostprocessing_drive_KEEN.f 
#YG                                  PROPERTIES COMPILE_FLAGS -qfixed)
#YG   ENDIF(Fortran_COMPILER MATCHES IBM)
#YG
#YG   ADD_EXECUTABLE(VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f)
#YG   TARGET_LINK_LIBRARIES(VPpostprocessing_drive_KEEN dfftpack)

ENDIF()

