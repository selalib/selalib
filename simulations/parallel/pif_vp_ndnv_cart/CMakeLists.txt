set(    SIM_NAME pif_vp_ndnv_cart             )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

add_executable( ${EXECUTABLE} sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
target_link_libraries( ${EXECUTABLE} ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} )

install( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
install( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
if(BUILD_TESTING)
  file( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  set( ARGS ${CMAKE_CURRENT_BINARY_DIR}/quick_test.nml )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
endif(BUILD_TESTING)
