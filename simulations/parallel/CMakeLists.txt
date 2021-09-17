# ##############################################################################
# Parallel simulations
# ##############################################################################
if(HDF5_ENABLED
   AND HDF5_IS_PARALLEL
   AND HDF5_PARALLEL_ENABLED)

  # Variables passed to simulations' CMakeFiles in subdirectories
  set(INSTALL_DIR "${INSTALL_DIR}/parallel")

  # Semi-Lagrangian simulations

  # Vlasov-Poisson 1Dx1D
  add_subdirectory(bsl_vp_1d1v_cart)
  add_subdirectory(bsl_vp_1d1v_cart_two_species)
  add_subdirectory(bsl_vp_1d1v_cart_multi_species)
  add_subdirectory(bsl_vp_1d1v_cart_no_split)

  # Vlasov-Ampere 1Dx1D
  add_subdirectory(bsl_va_1d1v_cart_berk_breizman)
  add_subdirectory(bsl_va_1d1v_cart_spectral)

  # Vlasov-Poisson 2Dx2D
  add_subdirectory(bsl_vp_2d2v_cart)
  add_subdirectory(bsl_vp_2d2v_cart_poisson_serial)
  if(NOT (FFT_LIB STREQUAL "SLLFFT"))
    add_subdirectory(bsl_vp_2d2v_cart_sparsegrid)
  endif()
  add_subdirectory(bsl_vp_2d2v_polar)
  add_subdirectory(bsl_vp_2d2v_curv)

  # Vlasov-Poisson 3Dx3D
  if(NOT (FFT_LIB STREQUAL "SLLFFT"))
    add_subdirectory(bsl_vp_3d3v_cart_sparsegrid)
  endif()

  # Drift/gyro-kinetic 3Dx1D
  add_subdirectory(bsl_dk_3d1v_polar)
  add_subdirectory(bsl_gk_3d1v_polar_one_mu)
  add_subdirectory(bsl_gk_3d1v_polar_multi_mu)
  add_subdirectory(bsl_dk_3d1v_polar_field_aligned)
  add_subdirectory(bsl_dk_3d1v_curv_field_aligned)

  # Eulerian simulations

  # Vlasov-Poisson 2Dx2D
  add_subdirectory(eul_vp_2d2v_cart_fv)

  # Particle simulations

  # Vlasov-Poisson 1Dx1D

  # Vlasov-Maxwell 1Dx2D
  add_subdirectory(pic_vm_1d2v_cart)
  add_subdirectory( pic_vm_1d2v_cart_multispecies )

  # Vlasov-Poisson 2Dx2D
  add_subdirectory(pic_vp_2d2v_cart)
  add_subdirectory(pic_vp_2d2v_cart_optim)
  add_subdirectory(pic_vp_2d2v_cart_optim_push)
  if(NOT (FFT_LIB STREQUAL "SLLFFT"))
    add_subdirectory(pic_vp_2d2v_cart_lbf)
  endif()

  # Vlasov-Poisson NDxND
  add_subdirectory(pif_vp_ndnv_cart)

  # Guiding-Center 2D
  add_subdirectory(pic_gc_2d0v_cart_optim)
  
  ## Vlasov-Maxwell 2Dx3D
  add_subdirectory( pic_vm_2d3v_cart )
  
  ## Vlasov-Maxwell 3Dx3D
  add_subdirectory( pic_vm_3d3v_cart )
  add_subdirectory( pic_vm_3d3v_cart_multispecies )
    


endif()
