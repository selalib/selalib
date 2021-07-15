set(    SIM_NAME pic_vm_1d2v_cart             )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(   TEST_NAME_HS "${TEST_NAME}_hs"            )
set(   TEST_NAME_HSSA "${TEST_NAME}_hs_stronga" )
set(   TEST_NAME_CEF   "${TEST_NAME}_cef"     )
set(   TEST_NAME_ECSIM "${TEST_NAME}_ecsim"   )
set(   TEST_NAME_disgradE   "${TEST_NAME}_disgradE"     )
set(   TEST_NAME_DISGRAD   "${TEST_NAME}_disgrad"     )
set(   TEST_NAME_DISGRAD_SA   "${TEST_NAME}_disgrad_stronga"     )
set(   TEST_NAME_DISGRADSUB   "${TEST_NAME}_disgradsub"     )
set(   TEST_NAME_hs_SMOOTH   "${TEST_NAME}_hs_smooth"     )
set(   TEST_NAME_disgradE_SMOOTH   "${TEST_NAME_disgradE}_smooth"     )
set(   TEST_NAME_DISGRAD_SMOOTH   "${TEST_NAME}_disgrad_smooth"     )
set(   TEST_NAME_DISGRADSUB_SMOOTH   "${TEST_NAME}_disgradsub_smooth"     )
set(   TEST_NAME_TRAFO "${TEST_NAME}_trafo"   )
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
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_hs.nml )
  add_mpi_test( ${TEST_NAME_HS} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME_HS} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )
    
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_hs_stronga.nml )
  add_mpi_test( ${TEST_NAME_HSSA} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME_HSSA} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )

  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_cef.nml )
  add_mpi_test( ${TEST_NAME_CEF} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_CEF} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )
  
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_ecsim.nml )
  add_mpi_test( ${TEST_NAME_ECSIM} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_ECSIM} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )

    set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_disgradE.nml )
  add_mpi_test( ${TEST_NAME_disgradE} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_disgradE} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )

  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_disgrad.nml )
  add_mpi_test( ${TEST_NAME_DISGRAD} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_DISGRAD} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )
    
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_disgrad_stronga.nml )
  add_mpi_test( ${TEST_NAME_DISGRAD_SA} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME_DISGRAD_SA} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )
  
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_disgradsub.nml )
  add_mpi_test( ${TEST_NAME_DISGRADSUB} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_DISGRADSUB} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 50 )

  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_hs_smooth.nml )
  add_mpi_test( ${TEST_NAME_hs_SMOOTH} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_hs_SMOOTH} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 20 )

  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_disgradE_smooth.nml )
  add_mpi_test( ${TEST_NAME_disgradE_SMOOTH} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_disgradE_SMOOTH} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 20 )

  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_disgrad_smooth.nml )
  add_mpi_test( ${TEST_NAME_DISGRAD_SMOOTH} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_DISGRAD_SMOOTH} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 20 )
  
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_disgradsub_smooth.nml )
  add_mpi_test( ${TEST_NAME_DISGRADSUB_SMOOTH} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_DISGRADSUB_SMOOTH} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 50 )

  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_trafo.nml )
  add_mpi_test( ${TEST_NAME_TRAFO} ${EXECUTABLE} ${PROCS} ${NML_FILE}  )
  set_tests_properties( ${TEST_NAME_TRAFO} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )

endif(BUILD_TESTING)
