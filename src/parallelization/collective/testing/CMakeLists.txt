set(ARGS " ")

add_executable(test_collective test_collective.F90)
target_link_libraries(test_collective sll_collective sll_timer)
add_mpi_test(collective test_collective ${PROCS} ${ARGS})
set_tests_properties(collective PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS")

# Remove this test because it raises errors on Fedora with openmpi
# ADD_EXECUTABLE( test_hello_mpi hello_mpi.F90 ) TARGET_LINK_LIBRARIES(
# test_hello_mpi sll_collective sll_timer ) SET(PROCS 4) ADD_MPI_TEST( hello_mpi
# test_hello_mpi ${PROCS} ${ARGS} ) SET_TESTS_PROPERTIES( hello_mpi PROPERTIES
# FAIL_REGULAR_EXPRESSION "NOT PASS" )
