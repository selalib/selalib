set(    SIM_NAME pif_vp_ndnv_cart             )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} 
  sll_sim_parallel
  sll_pif
  sll_timer
  sll_memory
  sll_utilities
  sll_constants
  sll_sobol
  sll_prob
  sll_pic_visu
  sll_pic_visu_parallel
  sll_collective
  sll_moment_matching
  sll_pic_utilities
  sll_descriptors
  sll_splitting_methods)


INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( ARGS ${CMAKE_CURRENT_BINARY_DIR}/quick_test.nml )
  ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)

