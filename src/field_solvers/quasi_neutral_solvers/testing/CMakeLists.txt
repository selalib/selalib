add_executable( test_qn_solver_2d_polar
  test_qn_solver_2d_polar.F90
  m_test_case_qn_2d_base.F90
  m_test_case_qn_2d_dirichlet.F90
  m_test_case_qn_2d_neumann_mode0.F90 )

target_link_libraries( test_qn_solver_2d_polar sll_quasi_neutral_solvers )

add_test( NAME qn_solver_2d_polar COMMAND test_qn_solver_2d_polar )
set_tests_properties( qn_solver_2d_polar
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 10 )
