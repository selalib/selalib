set(   TEST_NAME 2d_gc_curvilinear             )
set(    SIM_NAME 2d_guiding_center_curvilinear )
set( INPUT_FILES gcsim2d_curvilinear_input.nml 
                 gcsim2d_curvilinear_input_cart.nml
                 gcsim2d_curvilinear_input_polar.nml
                 gcsim2d_curvilinear_input_DSG.nml   )
set(   POST_PROC gcsim2d_cartesian_ref.dat
                 gcsim2d_curvilinear_ref.dat
                 gcsim2d_curvilinear.gnu 
                 gcsim2d_curvilinear_cart.gnu
                 gcsim2d_curvilinear_polar.gnu)

set(  EXECUTABLE "test_${TEST_NAME}" )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_sequential ${SIM_LIBS} )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

