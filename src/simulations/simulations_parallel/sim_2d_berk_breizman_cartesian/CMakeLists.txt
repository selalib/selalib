IF(FFTW_FOUND AND FFTW_ENABLED)

  FILE( GLOB INPUT_FILES "*.nml"         )
  FILE( GLOB   POST_PROC "*.dat" "*.gnu" )
  
  SET(   TEST_NAME 1d1d_vlasov_ampere )
  SET(    SIM_NAME 1d1d_vlasov_ampere )
  SET(  EXECUTABLE "test_${TEST_NAME}"           )
  SET( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )
  
  ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                                unit_test_${TEST_NAME}.F90 )
    
  TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel )
  
  INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
  INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )
  
ENDIF(FFTW_FOUND AND FFTW_ENABLED)
  
