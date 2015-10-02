FILE( GLOB INPUT_FILES "*.nml"         )
FILE( GLOB   POST_PROC "*.dat" "*.gnu" )

SET(   TEST_NAME 2d_vp_cartesian               )
SET(    SIM_NAME 2d_vlasov_poisson_cartesian   )
SET(  EXECUTABLE "test_${TEST_NAME}"           )
SET( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
  
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET(  ARGS ${CMAKE_CURRENT_BINARY_DIR}/vpsim2d_cartesian_input )
  ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)

#==============================================================================
# Postprocessing code for KEEN
#==============================================================================
IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES( VPpostprocessing_drive_KEEN.f 
                                PROPERTIES COMPILE_FLAGS -qfixed )
ENDIF(Fortran_COMPILER MATCHES IBM)

ADD_EXECUTABLE( VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f )
TARGET_LINK_LIBRARIES( VPpostprocessing_drive_KEEN dfftpack )

INSTALL( TARGETS VPpostprocessing_drive_KEEN DESTINATION ${INSTALL_DIR} )
INSTALL( FILES vp_plots.py DESTINATION ${INSTALL_DIR} )

IF(FFTW_FOUND AND FFTW_ENABLED)

  #PN Vlasov Ampere version of Michel program
  SET( TEST_NAME   2d_va_cartesian            )
  SET( SIM_NAME    2d_vlasov_ampere_cartesian )
  SET( EXECUTABLE "test_${TEST_NAME}"         )
  
  ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                                unit_test_${TEST_NAME}.F90 )
    
  TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel )
  
  INSTALL( TARGETS ${EXECUTABLE} DESTINATION ${INSTALL_DIR} )
  
  IF(BUILD_TESTING)
    FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
    SET(  ARGS ${CMAKE_CURRENT_BINARY_DIR}/vasim2d_cartesian_input )
    ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
    SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  ENDIF(BUILD_TESTING)
  
  #ADD_EXECUTABLE( vm2d_spectral sll_vlasov2d_base.F90
  #                              sll_vlasov2d_spectral.F90
  #                              init_functions.F90
  #                              vm2d_spectral.F90)
  #  
  #TARGET_LINK_LIBRARIES( vm2d_spectral sll_ampere_solvers )

ENDIF(FFTW_FOUND AND FFTW_ENABLED)
  
#================================================================
# add two species simulations
# for the moment, only a copy of one specy simulation

#comment for the moment as build fail on prototype-devel
SET(   TEST_NAME 2d_vp_cartesian_two_species               )
SET(    SIM_NAME 2d_vlasov_poisson_cartesian_two_species   )
SET(  EXECUTABLE "test_${TEST_NAME}"           )
SET( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
  
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel )

