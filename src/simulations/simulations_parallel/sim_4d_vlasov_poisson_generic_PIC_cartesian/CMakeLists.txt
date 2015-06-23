
set(   TEST_NAME 4d_vp_generic_pic_cartesian           )
set(    SIM_NAME 4d_vp_generic_pic_cartesian           )
set(  EXECUTABLE "test_${TEST_NAME}"           )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE}
  sll_simulations_parallel
  sll_working_precision
  sll_pic_interface
  sll_simple_pic_particles )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( ARGS ${CMAKE_CURRENT_BINARY_DIR}/params_generic_pic_4d.nml )
  ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
    TIMEOUT 200 )
ENDIF(BUILD_TESTING)

