set(   TEST_NAME 2d_gc_polar             )
set(    SIM_NAME 2d_guiding_center_polar )
set(  INPUT_FILE gcsim2d_polar_input.nml )
set(   POST_PROC gcsim2d_polar.gnu gcsim2d_polar_ref.dat )

set(  EXECUTABLE "test_${TEST_NAME}" )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_sequential ${SIM_LIBS} )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   ${INPUT_FILE} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Quick fix for ctest
FILE( COPY ${INPUT_FILE} DESTINATION ${CMAKE_BINARY_DIR} )

