set(    SIM_NAME pic_vm_3d3v_cart             )
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
  set(   TEST_NAME "sim_${SIM_NAME}_hs"            )
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_hs_testing.nml )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )

  # Test with alternative Hamiltonian splitting from crouseilles, einkemmer, faou
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_cef_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_cef"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )

  # Test with alternative Hamiltonian splitting from crouseilles, einkemmer, faou
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_cl_cef_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_cl_cef"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )

  # Test with average vector field
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_disgradE_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_disgradE"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )

   # Test with average vector field and boundary conditions
   set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_cl_disgradE_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_cl_disgradE"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} ) 
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )

  # Test with modified average vector field (truncated particle mass)
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_disgradE_trunc_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_disgradE_trunc"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )

  # Test with discrete gradient
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_disgradEC_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_disgradEC"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 150 )

  # Test with coordinate transformation and energy conserving splitting
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_disgradE_trafo_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_disgradE_trafo"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 500 )

   # Test with coordinate transformation and energy conserving splitting
   set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_cl_disgradE_trafo_testing.nml )
  set(   TEST_NAME "sim_${SIM_NAME}_cl_disgradE_trafo"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} ) 
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 500 )

  # Test with coordinate transformation and GEMPIC Hamiltonian splitting
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_hs_trafo_testing.nml )

  set(   TEST_NAME "sim_${SIM_NAME}_hs_trafo"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 500 )

  # Test with coordinate transformation and GEMPIC Hamiltonian splitting and boundary conditions
  #set( REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/reffile_pic_vm_3d3v_cart_cl_hs_trafo.dat )
  #file( COPY ${REFERENCE_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_cl_hs_trafo_testing.nml )

  set(   TEST_NAME "sim_${SIM_NAME}_cl_hs_trafo"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE})# ${REFERENCE_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 500 )

  # Test with coordinate transformation and cef Hamiltonian splitting
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_cef_trafo_testing.nml )
  
  set(   TEST_NAME "sim_${SIM_NAME}_cef_trafo"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 500 )

   # Test with coordinate transformation and cef Hamiltonian splitting and boundary conditions
   set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_cl_cef_trafo_testing.nml )
  
  set(   TEST_NAME "sim_${SIM_NAME}_cl_cef_trafo"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 500 )

  # Test with coordinate transformation and discrete gradient discretisation
  set( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_3d3v_cart_disgradEC_trafo_testing.nml )

  set(   TEST_NAME "sim_${SIM_NAME}_disgradEC_trafo"            )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 500 )
  
endif(BUILD_TESTING)

