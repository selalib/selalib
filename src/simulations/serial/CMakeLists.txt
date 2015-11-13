##############################################################################
###      Serial simulations
##############################################################################

ADD_LIBRARY( sll_sim_serial STATIC ../simulation_base_class.F90 )

TARGET_LINK_LIBRARIES( sll_sim_serial
                       sll_interpolators
                       sll_fields
                       sll_poisson_solvers
                       sll_reduction
                       sll_splitting_methods
                       sll_common_array_initializers
                       sll_gyroaverage
                     )

SET( INSTALL_DIR "${INSTALL_DIR}/sim_serial" )


IF(HDF5_ENABLED)
  
  # Semi Lagrangian method
  
  ## advection
  ADD_SUBDIRECTORY( bsl_ad_2d0v_cart                          )
  ADD_SUBDIRECTORY( bsl_ad_2d0v_polar_rot                     )
  ADD_SUBDIRECTORY( bsl_ad_2d0v_curv                          )
  ADD_SUBDIRECTORY( bsl_ad_2d0v_hex_rot                       )
  ADD_SUBDIRECTORY( fsl_ad_2d0v_polar                         )
  ADD_SUBDIRECTORY( fsl_ad_2d0v_curv                          )

  ## Vlasov-Poisson 1Dx1D
  ADD_SUBDIRECTORY( bsl_vp_1d1v_cart_serial                   )
  ADD_SUBDIRECTORY( bsl_vp_1d1v_cart_deltaf                   )
  ADD_SUBDIRECTORY( bsl_vp_1d1v_cart_micro_macro              )

  ## Vlasov-Poisson 2Dx2D
  ADD_SUBDIRECTORY( bsl_vp_2d2v_cart_fullsparsegrid           )

  ## Guiding-center 2D
  ADD_SUBDIRECTORY( bsl_gc_2d0v_cart                          )
  ADD_SUBDIRECTORY( bsl_gc_2d0v_polar                         )
  ADD_SUBDIRECTORY( bsl_gc_2d0v_polar_one_mu                  )
  ADD_SUBDIRECTORY( bsl_gc_2d0v_hex                           )
  ADD_SUBDIRECTORY( bsl_gc_2d0v_hex_splines                   )
  ADD_SUBDIRECTORY( bsl_gc_2d0v_hex_hermite                   )

  # Eulerian method

  ## Vlasov-Poisson 1Dx1D
  ADD_SUBDIRECTORY( eul_vp_1d1v_cart_fe_ft                    )

ENDIF(HDF5_ENABLED)

ADD_SUBDIRECTORY(bsl_examples)
ADD_SUBDIRECTORY(fsl_examples)
