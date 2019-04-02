# Save .mod files locally (this does not pollute the main "modules" directory)
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

add_executable( test_qn_solver_2d_polar
  test_qn_solver_2d_polar.F90
  m_test_qn_solver_2d_polar_base.F90
  m_test_qn_solver_2d_polar_annulus_dirichlet.F90
  m_test_qn_solver_2d_polar_annulus_neumann_mode0.F90
  m_test_qn_solver_2d_polar_disk_dirichlet.F90 )

target_link_libraries( test_qn_solver_2d_polar sll_quasi_neutral_solvers )

add_test( NAME qn_solver_2d_polar COMMAND test_qn_solver_2d_polar )
set_tests_properties( qn_solver_2d_polar
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 10 )

# Test FE quasi-neutrality solver based on smooth polar splines
add_executable( test_qn_solver_2d_fem_sps_stencil test_qn_solver_2d_fem_sps_stencil.F90 )
target_link_libraries( test_qn_solver_2d_fem_sps_stencil sll_quasi_neutral_solvers )
