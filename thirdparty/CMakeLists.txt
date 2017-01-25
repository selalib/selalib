#############################################################################
###      Third party simulations
#############################################################################

if(BUILD_THIRDPARTY)
  if( HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED )
    add_subdirectory( pic_vp_2d2v_cart_tau )
    add_subdirectory( bsl_vm_2d2v_cart_spectral/src )
  endif()
endif(BUILD_THIRDPARTY)

