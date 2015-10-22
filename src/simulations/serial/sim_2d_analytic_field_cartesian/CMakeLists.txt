set(   TEST_NAME 2d_af_cartesian             )
set(    SIM_NAME 2d_analytic_field_cartesian )
set(  INPUT_FILE afsim2d_cartesian_sdf.nml   )
set(   POST_PROC afsim2d_cartesian_sdf.gnu   )

set(  EXECUTABLE "test_${TEST_NAME}" )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_sim_serial ${SIM_LIBS} )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   ${INPUT_FILE} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )
