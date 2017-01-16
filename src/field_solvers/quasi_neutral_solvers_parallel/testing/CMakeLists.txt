add_executable( test_qn_solver_2d_polar_par
  test_qn_solver_2d_polar_par.F90
  m_test_case_qn_par_2d_base.F90
  m_test_case_qn_par_2d_dirichlet.F90
  m_test_case_qn_par_2d_neumann_mode0.F90 )

target_link_libraries( test_qn_solver_2d_polar_par
  sll_quasi_neutral_solvers_parallel )

ADD_MPI_TEST( qn_solver_2d_polar_par__np=1 test_qn_solver_2d_polar_par 1 " " )
ADD_MPI_TEST( qn_solver_2d_polar_par__np=4 test_qn_solver_2d_polar_par 4 " " )
set_tests_properties(
  qn_solver_2d_polar_par__np=1
  qn_solver_2d_polar_par__np=4
  PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 10 )
