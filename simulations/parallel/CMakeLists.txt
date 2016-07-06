#############################################################################
###      Parallel simulations
#############################################################################
if( HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED )

    # Variables passed to simulations' CMakeFiles in subdirectories
    set( INSTALL_DIR "${INSTALL_DIR}/parallel" )

    # Semi-Lagrangian simulations
    
    ## Vlasov-Poisson 1Dx1D
    add_subdirectory( bsl_vp_1d1v_cart                         )
    add_subdirectory( bsl_vp_1d1v_cart_two_species             )
    add_subdirectory( bsl_vp_1d1v_cart_multi_species           )
    add_subdirectory( bsl_vp_1d1v_cart_no_split                )

    ## Vlasov-Ampere 1Dx1D
    add_subdirectory( bsl_va_1d1v_cart_berk_breizman           )
    add_subdirectory( bsl_va_1d1v_cart_spectral                )

    ## Vlasov-Poisson 2Dx2D
    add_subdirectory( bsl_vp_2d2v_cart                         )
    add_subdirectory( bsl_vp_2d2v_cart_poisson_serial          )
    add_subdirectory( bsl_vp_2d2v_cart_sparsegrid              )
    add_subdirectory( bsl_vp_2d2v_cart_multipatch              )
    add_subdirectory( bsl_vp_2d2v_polar                        )
    add_subdirectory( bsl_vp_2d2v_curv                         )

    ## Vlasov-Poisson 3Dx3D
    add_subdirectory( bsl_vp_3d3v_cart_sparsegrid              )

    ## Guiding-Center 2D
    add_subdirectory( bsl_gc_2d0v_curv )

    ## Drift/gyro-kinetic 3Dx1D
    add_subdirectory( bsl_dk_3d1v_polar                        )
    add_subdirectory( bsl_gk_3d1v_polar_one_mu                 )
    add_subdirectory( bsl_gk_3d1v_polar_multi_mu               )
    if(FULL_FORTRAN2003)
      add_subdirectory( bsl_dk_3d1v_polar_field_aligned          )
    endif(FULL_FORTRAN2003)
    add_subdirectory( bsl_dk_3d1v_curv                         )
    add_subdirectory( bsl_dk_3d1v_curv_field_aligned           )


    # Eulerian simulations

    ## Vlasov-Poisson 2Dx2D
    if(FULL_FORTRAN2003)
      add_subdirectory( eul_vp_2d2v_cart_fv                      )
    endif(FULL_FORTRAN2003)

    # Particle simulations

    ## Vlasov-Poisson 1Dx1D

    ## Vlasov-Maxwell 1Dx2D
    add_subdirectory( pic_vm_1d2v_cart                         )

    ## Vlasov-Maxwell 2Dx2D
    ADD_SUBDIRECTORY( pic_vm_2d2v_cart                         )

    ## Vlasov-Poisson 2Dx2D
    add_subdirectory( pic_vp_2d2v_cart                         )
    add_subdirectory( pic_vp_2d2v_cart_optim                   )
    add_subdirectory( pic_vp_2d2v_cart_optim_push              )
    add_subdirectory( pic_vp_2d2v_cart_lbf                     )

    ## Vlasov-Poisson NDxND
    add_subdirectory( pif_vp_ndnv_cart                         )

    ## Guiding-Center 2D
    add_subdirectory( pic_gc_2d0v_cart_optim                   )
    
    
    #To delete/put somewhere else
    #as the simulation seems not to do what
    #is indicated
    #ADD_SUBDIRECTORY( sim_4d_drift_kinetic_cartesian          )
 
endif()


# here are some other examples of simulations
# they can serve as examples for selalib beginners
# you can add some new examples
# the interface is more free
# requirements of the simulations that are
# present in the developer guide are not
# mandatory here
add_subdirectory(bsl_examples)
