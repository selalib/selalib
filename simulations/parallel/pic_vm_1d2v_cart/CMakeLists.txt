if( FULL_FORTRAN2003 )

set(    SIM_NAME pic_vm_1d2v_cart             )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

add_executable( ${EXECUTABLE} sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
target_link_libraries( ${EXECUTABLE} ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} )

install( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
install( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

## Ctest
if(BUILD_TESTING)
  set( INPUT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/reffile_pic_vm_1d2v_cart.dat )
  file( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_testing.nml )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} ${REFERENCE_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )
endif(BUILD_TESTING)

endif( FULL_FORTRAN2003 )
