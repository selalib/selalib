set(   TEST_NAME general_pif                   )
set(    SIM_NAME general_vlasov_poisson_pif    )
set(  EXECUTABLE "test_${TEST_NAME}"            )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}"  )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE}        sll_simulations_parallel
					      sll_pif_fieldsolver
					      sll_timer
					      sll_memory
					      sll_utilities
					      sll_constants
					      sll_sobol
					      sll_prob
					      sll_visu_pic
					      sll_visu_pic_coll
					      sll_collective
					      sll_moment_matching
					      sll_pic_utilities
					      sll_wedge_product_generaldim
					      sll_descriptors
					      sll_splitting_methods
					      sll_general_vlasov_poisson_pif)

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( ARGS ${CMAKE_CURRENT_BINARY_DIR}/quick_test.nml )
  ADD_MPI_TEST( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)


					      
# Ctest
# IF(BUILD_TESTING)
#   ADD_EXECUTABLE( test_pic_accumulator unit_test.F90 )
#   TARGET_LINK_LIBRARIES( test_pic_accumulator
#     sll_pic_utilities
#     sll_pic_initializers
#     sll_pic_accumulators
#     sll_pic_particles
#     sll_meshes
#     sll_assert )
#   ADD_TEST( NAME pic_accumulator COMMAND test_pic_accumulator )
#   SET_TESTS_PROPERTIES( pic_accumulator PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
# ENDIF(BUILD_TESTING)


 
