set(   TEST_NAME 2d_gc_curvilinear             )
set(    SIM_NAME 2d_guiding_center_curvilinear )
set(  EXECUTABLE "test_${TEST_NAME}"           )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

