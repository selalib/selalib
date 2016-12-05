##############################################################################
###      Serial simulations
##############################################################################

set( INSTALL_DIR "${INSTALL_DIR}/serial" )

if(HDF5_ENABLED)

  # Semi Lagrangian method

  ## advection
  add_subdirectory( bsl_ad_2d0v_cart                          )
  add_subdirectory( bsl_ad_2d0v_polar_rot                     )
  add_subdirectory( bsl_ad_2d0v_curv                          )
  add_subdirectory( bsl_ad_2d0v_hex_rot                       )
  add_subdirectory( bsl_ad_2d0v_hex                           )
  add_subdirectory( fsl_ad_2d0v_polar                         )
  add_subdirectory( fsl_ad_2d0v_curv                          )

  ## Vlasov-Poisson 1Dx1D
  add_subdirectory( bsl_vp_1d1v_cart_serial                   )
  add_subdirectory( bsl_vp_1d1v_cart_deltaf                   )
  add_subdirectory( bsl_vp_1d1v_cart_micro_macro              )
  add_subdirectory( fsl_vp_1d1v_cart_tau                      )

  ## Vlasov-Poisson 2Dx2D
  add_subdirectory( bsl_vp_2d2v_cart_fullsparsegrid           )

  ## Guiding-center 2D
  add_subdirectory( bsl_gc_2d0v_cart                          )
  add_subdirectory( bsl_gc_2d0v_polar                         )
  add_subdirectory( bsl_gc_2d0v_polar_one_mu                  )
  add_subdirectory( bsl_gc_2d0v_hex                           )
  add_subdirectory( bsl_gc_2d0v_hex_splines                   )
  add_subdirectory( bsl_gc_2d0v_hex_hermite                   )

  # Eulerian method

  ## Vlasov-Poisson 1Dx1D
  add_subdirectory( eul_vp_1d1v_cart_fe_ft                    )

endif(HDF5_ENABLED)

# Examples that do not extend 'sll_c_simulation_base_class'
add_subdirectory( bsl_examples )
add_subdirectory( fsl_examples )
