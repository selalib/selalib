add_library(
  sll_linear_solvers STATIC
  sll_m_tridiagonal.F90
  sll_m_penta_diagonal.F90
  sll_m_choleski.F90
  sol.f
  muls.f
  sll_m_pivotbande.F90
  sll_m_linear_operator_base.F90
  sll_m_linear_operator_matrix_dense_to_dense.F90
  sll_m_linear_operator_matrix_dense_to_stencil.F90
  sll_m_linear_operator_matrix_dense_to_stencil_new.F90
  sll_m_linear_operator_matrix_stencil_to_dense.F90
  sll_m_linear_operator_matrix_stencil_to_dense_new.F90
  sll_m_linear_operator_matrix_stencil_to_stencil.F90
  sll_m_linear_operator_matrix_c1_block.F90
  sll_m_linear_operator_matrix_c1_block_new.F90
  sll_m_conjugate_gradient.F90
  sll_m_linear_operator_abstract.F90
  sll_m_matrix_abstract.F90
  sll_m_matrix_csr.F90
  sll_m_linear_operator_kron.F90
  sll_m_linear_operator_block.F90
  sll_m_linear_operator_penalized.F90
  sll_m_linear_solver_abstract.F90
  sll_m_linear_solver_iter_abstract.F90
  sll_m_linear_solver_kron.F90
  sll_m_linear_solver_block.F90
  sll_m_linear_solver_cg.F90
  sll_m_linear_solver_mgmres.F90
  sll_m_uzawa_iterator.F90
  sll_m_linear_operator_base.F90)

if(Fortran_COMPILER MATCHES IBM)
  set_source_files_properties(sol.f muls.f PROPERTIES COMPILE_FLAGS '-qfixed')
endif()

target_link_libraries(sll_linear_solvers sll_vector_space sll_errors sll_memory
                      sll_assert sll_sparsekit)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
