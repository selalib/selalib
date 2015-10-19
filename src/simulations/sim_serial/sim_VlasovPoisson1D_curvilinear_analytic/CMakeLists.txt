set(   TEST_NAME VP1D_curvilinear_analytic )
set(    SIM_NAME VP1D_curvilinear_analytic )

set(  EXECUTABLE "test_${TEST_NAME}" )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

IF(MPI_ENABLED)

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_sim_serial ${SIM_LIBS} )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   ${INPUT_FILE} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

ENDIF(MPI_ENABLED)
