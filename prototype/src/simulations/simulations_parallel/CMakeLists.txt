#############################################################################
###      Parallel simulations
#############################################################################

IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

   ADD_LIBRARY( sll_simulations_parallel STATIC ../simulation_base_class.F90)

   #Add new dependencies concerning simulation_MYMETHOD.F90
   TARGET_LINK_LIBRARIES(sll_simulations_parallel
                         sll_distribution_function_multipatch
                         sll_general_coordinate_elliptic_solvers_multipatch
                         sll_scalar_field_multipatch
                         sll_parallel_array_initializers
                         sll_advection
                         sll_parallel_utilities
                         sll_poisson_solvers
                         sll_time_splitting
                         sll_file_io_parallel
                         sll_poisson_parallel
                         sll_scalar_field_1d_alternative
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
                         sll_scalar_field_2d_alternative
                         sll_pic_particles
                         sll_random_deviate_generators
                         sll_field_accumulators
                         sll_scalar_field
                         sll_interpolators
                         sll_point_to_point_communications
                         sll_splines
                         sll_remap
                         sll_mesh_calculus
                         sll_integration
                         sll_constants
                         sll_coordinate_transformations
                         sll_logical_meshes
                         sll_utilities
                         sll_timer
                         sll_memory
                         sll_assert 	
                         sll_working_precision
                        )

    # Variables passed to simulations' CMakeFiles
    SET( INSTALL_DIR "${INSTALL_DIR}/simulations_parallel" )
    SET( SIM_LIBS sll_simulations_parallel )

    # Simulations as new subdirectories
    ADD_SUBDIRECTORY( sim_2d_guiding_center_curvilinear  )
    ADD_SUBDIRECTORY( sim_2d_vlasov_poisson_cartesian    )
    ADD_SUBDIRECTORY( sim_2d_vlasov_poisson_no_splitting )

#   MACRO(ADD_SIMULATION TEST_NAME SIM_NAME) 
#      ADD_EXECUTABLE(test_${TEST_NAME} simulation_${SIM_NAME}.F90
#                                       unit_test_${TEST_NAME}.F90)
#      TARGET_LINK_LIBRARIES(test_${TEST_NAME} sll_simulations)
#   ENDMACRO(ADD_SIMULATION)
#
#   #Define the executables/tests here.
#   #Please add your new test_MYMETHOD executable, built on unit_test_MYMETHOD.F90
#   #and simulation_MYSIM.F90
#   #ADD_SIMULATION(MYMETHOD MYSIM)
#
#   ADD_SIMULATION(2d_vp_cartesian   2d_vlasov_poisson_cartesian)
#   ADD_SIMULATION(2d_gc_curvilinear 2d_guiding_center_curvilinear)
#   ADD_SIMULATION(4d                4d_vlasov_poisson_cartesian)
#   ADD_SIMULATION(4d_vp_cartesian   4d_vlasov_parallel_poisson_sequential_cartesian)
#   ADD_SIMULATION(4d_vp_general     4d_vlasov_poisson_general)
#   ADD_SIMULATION(4d_qns_general    4d_qns_general)
#   ADD_SIMULATION(4d_qns_general_multipatch    4d_qns_general_multipatch)
#   ADD_SIMULATION(4d_vp_polar       4d_vlasov_poisson_polar)
#
#   ADD_SIMULATION(4d_DK_hybrid              4d_DK_hybrid)
#   ADD_SIMULATION(4d_dk_polar               4d_drift_kinetic_polar)
#   ADD_SIMULATION(4d_dk_polar_one_mu        4d_drift_kinetic_polar_one_mu)
#   ADD_SIMULATION(4d_dk_polar_multi_mu      4d_drift_kinetic_polar_multi_mu)
#   ADD_SIMULATION(4d_dk_field_aligned_polar 4d_drift_kinetic_field_aligned_polar)
#   ADD_SIMULATION(4d_dk_cartesian           4d_drift_kinetic_cartesian_finite_volume)
#   ADD_SIMULATION(6d                        6d_vlasov_poisson_cartesian)
#
#   ADD_SIMULATION(4d_vp_pic_cartesian       4d_vp_pic_cartesian)
#
##   ADD_SIMULATION(4d_mixed_qns_cartesian 4d_mixed_qns_cartesian)
#
#   # Add postprocessing code for KEEN
#   IF(Fortran_COMPILER MATCHES IBM)
#      SET_SOURCE_FILES_PROPERTIES(VPpostprocessing_drive_KEEN.f 
#                                  PROPERTIES COMPILE_FLAGS -qfixed)
#   ENDIF(Fortran_COMPILER MATCHES IBM)
#
#   ADD_EXECUTABLE(VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f)
#   TARGET_LINK_LIBRARIES(VPpostprocessing_drive_KEEN dfftpack)

ENDIF()

