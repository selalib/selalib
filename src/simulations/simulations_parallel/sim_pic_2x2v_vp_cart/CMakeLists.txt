set ( TEST_NAME pic_2x2v_vp_cart )
set ( SIM_NAME pic_2x2v_vp_cart )
set ( EXECUTABLE "test_${TEST_NAME}")
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )

ADD_EXECUTABLE( ${EXECUTABLE} sim_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE}
  sll_simulations_parallel
  sll_working_precision
  sll_pic_interface
  sll_particle_initializers
  sll_particle_types
  sll_kernel_smoothers
  sll_particle_pushers
  sll_poisson_solvers)

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( ARGS ${CMAKE_CURRENT_BINARY_DIR}/params_pic_vp_2x2v.nml )
  ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 200 )
ENDIF(BUILD_TESTING)
