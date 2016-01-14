#############################################################################
###      Parallel simulations
#############################################################################
IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)

  ADD_LIBRARY( sll_sim_parallel STATIC ${SIM_BASE_CLASS})

   # Add new dependencies concerning SIMNAME/sim_SIMNAME.F90 
   # and or SIMNAME/sll_m_sim_SIMNAME.F90
   # refer to the developer guide for choosing the right SIMNAME

   TARGET_LINK_LIBRARIES(sll_sim_parallel
                         sll_collective
                         sll_parallel_array_utilities
                         sll_advection
                         sll_parallel_utilities
                         sll_poisson_solvers
                         sll_splitting_methods
                         sll_file_io_parallel
                         sll_poisson_parallel
                         sll_fcisl
                         sll_quasi_neutral_solvers_parallel
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
                         sll_quadrature
                         sll_coordinate_transformations
                         sll_meshes
                         sll_utilities
                         sll_timer
                         sll_memory
                         sll_assert
                         )

    # Variables passed to simulations' CMakeFiles in subdirectories
    SET( INSTALL_DIR "${INSTALL_DIR}/parallel" )

    # Semi-Lagrangian simulations
    
    ## Vlasov-Poisson 1Dx1D
    ADD_SUBDIRECTORY( bsl_vp_1d1v_cart                         )
    ADD_SUBDIRECTORY( bsl_vp_1d1v_cart_two_species             )
    ADD_SUBDIRECTORY( bsl_vp_1d1v_cart_multi_species           )
    ADD_SUBDIRECTORY( bsl_vp_1d1v_cart_no_split                )

    ## Vlasov-Ampere 1Dx1D
    ADD_SUBDIRECTORY( bsl_va_1d1v_cart_berk_breizman           )
    ADD_SUBDIRECTORY( bsl_va_1d1v_cart_spectral                )

    ## Vlasov-Poisson 2Dx2D
    ADD_SUBDIRECTORY( bsl_vp_2d2v_cart                         )
    ADD_SUBDIRECTORY( bsl_vp_2d2v_cart_poisson_serial          )
    ADD_SUBDIRECTORY( bsl_vp_2d2v_cart_sparsegrid              )
    ADD_SUBDIRECTORY( bsl_vp_2d2v_cart_multipatch              )
    ADD_SUBDIRECTORY( bsl_vp_2d2v_polar                        )
    ADD_SUBDIRECTORY( bsl_vp_2d2v_curv                         )

    ## Vlasov-Poisson 3Dx3D
    ADD_SUBDIRECTORY( bsl_vp_3d3v_cart                         )
    ADD_SUBDIRECTORY( bsl_vp_3d3v_cart_sparsegrid              )

    ## Guiding-Center 2D
    ADD_SUBDIRECTORY( bsl_gc_2d0v_curv )

    ## Drift/gyro-kinetic 3Dx1D
    ADD_SUBDIRECTORY( bsl_dk_3d1v_polar                        )
    ADD_SUBDIRECTORY( bsl_gk_3d1v_polar_one_mu                 )
    ADD_SUBDIRECTORY( bsl_gk_3d1v_polar_multi_mu               )
    IF(FULL_FORTRAN2003)
      ADD_SUBDIRECTORY( bsl_dk_3d1v_polar_field_aligned          )
    ENDIF(FULL_FORTRAN2003)
    ADD_SUBDIRECTORY( bsl_dk_3d1v_curv                         )
    ADD_SUBDIRECTORY( bsl_dk_3d1v_curv_field_aligned           )


    # Eulerian simulations

    ## Vlasov-Poisson 2Dx2D
    IF(FULL_FORTRAN2003)
    ADD_SUBDIRECTORY( eul_vp_2d2v_cart_fv                      )
    ENDIF(FULL_FORTRAN2003)

    # Particle simulations

    ## Vlasov-Poisson 1Dx1D
    ADD_SUBDIRECTORY( pic_vp_1d1v_cart                         )

    ## Vlasov-Maxwell 1Dx2D
    ADD_SUBDIRECTORY( pic_vm_1d2v_cart                         )

    ## Vlasov-Poisson 2Dx2D
    ADD_SUBDIRECTORY( pic_vp_2d2v_cart                         )
    ADD_SUBDIRECTORY( pic_vp_2d2v_cart_optim                   )
    ADD_SUBDIRECTORY( pic_vp_2d2v_cart_optim_push              )
    ADD_SUBDIRECTORY( pic_vp_2d2v_cart_remapped                )

    ## Vlasov-Poisson NDxND
    ADD_SUBDIRECTORY( pif_vp_ndnv_cart                         )

    ## Guiding-Center 2D
    ADD_SUBDIRECTORY( pic_gc_2d0v_cart_optim                   )
    
    
    #To delete/put somewhere else
    #as the simulation seems not to do what
    #is indicated
    #ADD_SUBDIRECTORY( sim_4d_drift_kinetic_cartesian          )
 
ENDIF()

IF(MPI_ENABLED)
   # here are some other examples of simulations
   # they can serve as examples for selalib beginners
   # you can add some new examples
   # the interface is more free
   # requirements of the simulations that are
   # present in the developer guide are not
   # mandatory here
   ADD_SUBDIRECTORY(bsl_examples)
ENDIF(MPI_ENABLED)
