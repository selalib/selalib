set(    SIM_NAME pic_vm_1d2v_cart             )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )


file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE}
  sll_sim_parallel
  sll_pic_interface
  sll_particle_types
  sll_kernel_smoothers
  sll_particle_pushers
  sll_maxwell_solvers
  sll_xdmf_io)

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# For the moment the simulation is not prepared to be a ctest
## Ctest
IF(BUILD_TESTING)
  SET( INPUT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/reffile_pic_vm_1d2v_cart.dat )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/params_pic_vm_1d2v_cart_testing.nml )
  ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${NML_FILE} ${REFERENCE_FILE} )
  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 5 )
ENDIF(BUILD_TESTING)
