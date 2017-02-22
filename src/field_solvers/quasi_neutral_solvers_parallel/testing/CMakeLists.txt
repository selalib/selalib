# Save .mod files locally (this does not pollute the main "modules" directory)
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

add_executable( test_qn_solver_2d_polar_par
  test_qn_solver_2d_polar_par.F90
  m_test_qn_solver_2d_polar_base.F90
  m_test_qn_solver_2d_polar_dirichlet.F90
  m_test_qn_solver_2d_polar_neumann_mode0.F90 )

target_link_libraries( test_qn_solver_2d_polar_par
  sll_quasi_neutral_solvers_parallel )

ADD_MPI_TEST( qn_solver_2d_polar_par_np1 test_qn_solver_2d_polar_par 1 " " )
ADD_MPI_TEST( qn_solver_2d_polar_par_np4 test_qn_solver_2d_polar_par 4 " " )
set_tests_properties(
  qn_solver_2d_polar_par_np1
  qn_solver_2d_polar_par_np4
  PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 10 )
