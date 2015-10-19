set(   TEST_NAME pic1d1v_vp_periodic           )
set(    SIM_NAME pic1d1v_vp_periodic           )
set(  EXECUTABLE "test_${TEST_NAME}"           )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" "*.py" )

ADD_EXECUTABLE( ${EXECUTABLE} sim_${SIM_NAME}.F90 unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_sim_parallel sll_pic1d )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )
