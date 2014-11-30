ADD_LIBRARY( sll_collective STATIC sll_collective.F90 )
TARGET_LINK_LIBRARIES( sll_collective sll_utilities )

# Ctests
IF(BUILD_TESTING)
  SET( ARGS " " )

  ADD_EXECUTABLE( test_collective unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_collective sll_collective sll_timer )
  ADD_MPI_TEST( collective test_collective ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( collective PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS" )

  #Remove this test because it raises errors on Fedora with openmpi
  #ADD_EXECUTABLE( test_hello_mpi hello_mpi.F90 )
  #TARGET_LINK_LIBRARIES( test_hello_mpi sll_collective sll_timer )
  #SET(PROCS 4)
  #ADD_MPI_TEST( hello_mpi test_hello_mpi ${PROCS} ${ARGS} )
  #SET_TESTS_PROPERTIES( hello_mpi PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS" )

ENDIF(BUILD_TESTING)
