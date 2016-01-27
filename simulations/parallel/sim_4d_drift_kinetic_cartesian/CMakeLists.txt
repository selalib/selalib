set(   TEST_NAME 4d_dk_cartesian                )
set(    SIM_NAME 4d_drift_kinetic_cartesian_finite_volume )
set(  EXECUTABLE "test_${TEST_NAME}"            )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}"  )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

add_executable( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
target_link_libraries( ${EXECUTABLE} ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} )

install( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
install( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
if(BUILD_TESTING)
  set( ARGS " " )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
endif(BUILD_TESTING)
