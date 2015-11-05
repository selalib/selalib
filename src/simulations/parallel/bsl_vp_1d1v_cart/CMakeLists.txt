set(    SIM_NAME bsl_vp_1d1v_cart             )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_sim_parallel )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )



FILE( GLOB INPUT_FILES "*.nml"         )
FILE( GLOB   POST_PROC "*.dat" "*.gnu" )




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


