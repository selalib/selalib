#############################################################################
###      Parallel simulations
#############################################################################
IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)

   ADD_LIBRARY( sll_simulations_parallel STATIC ../simulation_base_class.F90)

   #Add new dependencies concerning simulation_MYMETHOD.F90
   TARGET_LINK_LIBRARIES(sll_simulations_parallel
                         sll_collective
                         sll_parallel_array_utilities
                         sll_advection
                         sll_parallel_utilities
                         sll_poisson_solvers
                         sll_time_splitting
                         sll_file_io_parallel
                         sll_poisson_parallel
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

    # Variables passed to simulations' CMakeFiles in subdirectories
    SET( INSTALL_DIR "${INSTALL_DIR}/simulations_parallel" )

    # Simulations as new subdirectories
    ADD_SUBDIRECTORY( sim_2d_guiding_center_curvilinear        )
    ADD_SUBDIRECTORY( sim_2d_vlasov_poisson_cartesian          )
    ADD_SUBDIRECTORY( sim_2d_vlasov_poisson_no_splitting       )
    ADD_SUBDIRECTORY( sim_4d                                   )
    ADD_SUBDIRECTORY( sim_4d_vlasov_poisson_cartesian          )
    ADD_SUBDIRECTORY( sim_4d_vlasov_poisson_polar              )
    ADD_SUBDIRECTORY( sim_4d_vlasov_poisson_general            )
    ADD_SUBDIRECTORY( sim_4d_drift_kinetic_cartesian           )
    ADD_SUBDIRECTORY( sim_4d_drift_kinetic_polar               )
    ADD_SUBDIRECTORY( sim_4d_drift_kinetic_polar_one_mu        )
    ADD_SUBDIRECTORY( sim_4d_drift_kinetic_polar_multi_mu      )
    ADD_SUBDIRECTORY( sim_4d_drift_kinetic_field_aligned_polar )
    ADD_SUBDIRECTORY( sim_4d_drift_kinetic_hybrid              )
    ADD_SUBDIRECTORY( sim_4d_mixed_qns_cartesian               )
    ADD_SUBDIRECTORY( sim_4d_qns_general                       )
    ADD_SUBDIRECTORY( sim_6d                                   )

    ADD_SUBDIRECTORY( sim_sl_va_1d1v_cart_berk_breizman        )
    ADD_SUBDIRECTORY( sim_sl_va_1d1v_cart_spectral             )
    ADD_SUBDIRECTORY( sim_sl_vp_2d2v_cart_sparsegrid           )
    ADD_SUBDIRECTORY( sim_sl_vp_3d3v_cart_sparsegrid           )
    ADD_SUBDIRECTORY( sim_sl_vp_2d2v_cart_multipatch           )
    ADD_SUBDIRECTORY( sl_dk_3d1v_curv_field_aligned        )
 

    # Particle simulations (mostly PIC)
    ADD_SUBDIRECTORY( sim_2d_guiding_center_PIC_cartesian      )
    ADD_SUBDIRECTORY( sim_4d_vlasov_poisson_PIC_cartesian      )
    ADD_SUBDIRECTORY( sim_4d_vlasov_poisson_generic_PIC_cartesian )
    ADD_SUBDIRECTORY( sim_general_vlasov_poisson_pif           )
    ADD_SUBDIRECTORY( sim_pic1d1v_vlasov-poisson_periodic      )

    ADD_SUBDIRECTORY( sim_pic_vp_2d2v_cart                     )
    ADD_SUBDIRECTORY( sim_pic_vm_1d2v_cart                     )

ENDIF()

IF(MPI_ENABLED)
   ADD_SUBDIRECTORY(bsl_examples)
ENDIF(MPI_ENABLED)
