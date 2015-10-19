IF (FFTW_FOUND AND FFTW_ENABLED)
  set( SIM_NAME  sim_sl_vp_2d2v_cart_sparsegrid    )
  set( EXECUTABLE "${SIM_NAME}"            )
  set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}"  )
  
  file( GLOB INPUT_FILES "*.nml"         )
  
  ADD_EXECUTABLE( ${EXECUTABLE} sll_m_${SIM_NAME}.F90
    ${SIM_NAME}.F90 )
  TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_sim_parallel sll_sparse_grid )
  
  INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
  INSTALL( FILES                                     DESTINATION ${INSTALL_DIR} )
ENDIF()
