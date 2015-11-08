set( SIM_NAME    bsl_vp_1d1v_cart_multi_species )
set( TEST_NAME   "sim_${SIM_NAME}"              )
set( EXECUTABLE  "${TEST_NAME}"                 )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}"   )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.ref" "*.gnu" "*.pdf")

ADD_EXECUTABLE( ${EXECUTABLE} sll_m_species.F90 sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_sim_parallel )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )
