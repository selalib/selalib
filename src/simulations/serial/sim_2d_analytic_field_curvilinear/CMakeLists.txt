set(   TEST_NAME 2d_af_curvilinear             )
set(    SIM_NAME 2d_analytic_field_curvilinear )
set(  INPUT_FILE afsim2d_curvilinear_sdf.nml )

set(  EXECUTABLE "test_${TEST_NAME}" )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_sequential ${SIM_LIBS} )

INSTALL( TARGETS ${EXECUTABLE}         DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   ${INPUT_FILE} README  DESTINATION ${INSTALL_DIR} )

