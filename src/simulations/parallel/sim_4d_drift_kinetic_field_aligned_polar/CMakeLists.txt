set(   TEST_NAME 4d_dk_field_aligned_polar            )
set(    SIM_NAME 4d_drift_kinetic_field_aligned_polar )
set(  EXECUTABLE "test_${TEST_NAME}"                  )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}"        )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE}
  sll_sim_parallel
  sll_xdmf_io_parallel )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( ARGS ${CMAKE_CURRENT_BINARY_DIR}/dksim4d_field_aligned_polar.nml )
  ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES 
    PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 100 )
ENDIF(BUILD_TESTING)
