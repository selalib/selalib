# Sequential test
ADD_EXECUTABLE( test_qn_solver_2d test_qn_solver_2d.F90 )
TARGET_LINK_LIBRARIES( test_qn_solver_2d sll_quasi_neutral_solvers )
ADD_TEST( NAME qn_solver_2d COMMAND test_qn_solver_2d )
SET_TESTS_PROPERTIES( qn_solver_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  