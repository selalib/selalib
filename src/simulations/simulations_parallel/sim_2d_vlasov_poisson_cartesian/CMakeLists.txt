set(   TEST_NAME 2d_vp_cartesian               )
set(    SIM_NAME 2d_vlasov_poisson_cartesian   )
set(  EXECUTABLE "test_${TEST_NAME}"           )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )

TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel )

IF(Fortran_COMPILER MATCHES gfortran)
  ADD_EXECUTABLE( vm2d_spectral sll_vlasov2d_base.F90
                                sll_vlasov2d_spectral.F90
                                init_functions.F90
                                vm2d_spectral.F90)

  TARGET_LINK_LIBRARIES( vm2d_spectral sll_simulations_parallel )
ENDIF(Fortran_COMPILER MATCHES gfortran)

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( ARGS ${CMAKE_CURRENT_BINARY_DIR}/vpsim2d_cartesian_input )
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


