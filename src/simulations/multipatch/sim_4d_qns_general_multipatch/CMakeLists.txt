set(   TEST_NAME 4d_qns_general_multipatch      )
set(    SIM_NAME 4d_qns_general_multipatch      )
set(  EXECUTABLE "test_${TEST_NAME}"            )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}"  )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel 
                         sll_distribution_function_multipatch
                         sll_general_coordinate_elliptic_solvers_multipatch
                         sll_scalar_field_multipatch)

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

