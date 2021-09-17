set(    SIM_NAME pic_vm_2d3v_cart             )
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
  # Test with hs Hamiltonian splitting
  set( REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/reffile_pic_vm_2d3v_cart_hs.dat )
  file( COPY ${REFERENCE_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_2d3v_cart_hs_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_sympl"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} ${REFERENCE_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )
  
endif(BUILD_TESTING)
