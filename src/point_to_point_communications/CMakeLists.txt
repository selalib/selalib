ADD_LIBRARY( sll_point_to_point_communications STATIC
  sll_point_to_point_comms.F90 )

TARGET_LINK_LIBRARIES( sll_point_to_point_communications
  sll_collective )


# Ctests
IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_p2p_comms_1d unit_test_1d.F90 )
  TARGET_LINK_LIBRARIES( test_p2p_comms_1d
    sll_point_to_point_communications 
    sll_collective 
    sll_timer )

  ADD_EXECUTABLE( test_p2p_comms_2d unit_test_2d.F90 )
  TARGET_LINK_LIBRARIES( test_p2p_comms_2d
    sll_point_to_point_communications 
    sll_collective 
    sll_timer )

  SET( ARGS " " )
  ADD_MPI_TEST( point_to_point_comms_1d test_p2p_comms_1d ${PROCS} ${ARGS} )
  ADD_MPI_TEST( point_to_point_comms_2d test_p2p_comms_2d ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( point_to_point_comms_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( point_to_point_comms_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" ) 

ENDIF(BUILD_TESTING)
