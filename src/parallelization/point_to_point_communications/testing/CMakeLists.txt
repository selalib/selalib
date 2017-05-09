ADD_EXECUTABLE( test_point_to_point_comms_1d test_point_to_point_comms_1d.F90 )
TARGET_LINK_LIBRARIES( test_point_to_point_comms_1d
  sll_point_to_point_communications 
  sll_collective 
  sll_timer )

ADD_EXECUTABLE( test_point_to_point_comms_2d test_point_to_point_comms_2d.F90 )
TARGET_LINK_LIBRARIES( test_point_to_point_comms_2d
  sll_point_to_point_communications 
  sll_collective 
  sll_timer )

SET( ARGS " " )
ADD_MPI_TEST( point_to_point_comms_1d test_point_to_point_comms_1d ${PROCS} ${ARGS} )
ADD_MPI_TEST( point_to_point_comms_2d test_point_to_point_comms_2d "4" ${ARGS} )
SET_TESTS_PROPERTIES( point_to_point_comms_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
SET_TESTS_PROPERTIES( point_to_point_comms_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" ) 
