set(   TEST_NAME guiding_center_2D_generalized_coords_mudpack )
set(    SIM_NAME guiding_center_2D_generalized_coords_mudpack )

set(  EXECUTABLE "test_${TEST_NAME}" )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_sequential ${SIM_LIBS} )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   ${INPUT_FILE} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

