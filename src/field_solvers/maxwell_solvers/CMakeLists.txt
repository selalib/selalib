set(MAXWELL_SOLVERS_SRCS
  sll_m_maxwell_solvers_base.F90
  sll_m_maxwell_1d_base.F90
  sll_m_maxwell_1d_fem.F90
  sll_m_maxwell_1d_ps.F90
  sll_m_maxwell_1d_fem_sm.F90
  sll_m_maxwell_1d_trafo.F90
  sll_m_maxwell_2d_fdtd.F90
  sll_m_maxwell_2d_diga.F90
  sll_m_maxwell_2d_fem_fft.F90
  sll_m_maxwell_3d_base.F90
  sll_m_maxwell_3d_fem_fft.F90
  sll_m_maxwell_3d_fem.F90
  sll_m_maxwell_3d_trafo.F90
  sll_m_maxwell_clamped_1d_fem_sm.F90
  sll_m_maxwell_clamped_1d_trafo.F90
  sll_m_maxwell_clamped_3d_fem.F90
  sll_m_maxwell_clamped_3d_trafo.F90
  sll_m_preconditioner_fft.F90
  sll_m_preconditioner_singular.F90
  sll_m_linear_operator_poisson_1d.F90
  sll_m_linear_operator_poisson_3d.F90
  sll_m_linear_operator_poisson_clamped_1d.F90
  sll_m_linear_operator_poisson_clamped_3d.F90
  sll_m_linear_operator_schur_eb_1d.F90
  sll_m_linear_operator_schur_eb_3d.F90
  sll_m_linear_operator_schur_eb_cl_1d.F90
  sll_m_linear_operator_schur_eb_cl_3d.F90
  sll_m_linear_operator_maxwell_eb_schur.F90)

list(APPEND MAXWELL_SOLVERS_SRCS sll_m_maxwell_2d_pstd.F90
     sll_m_maxwell_3d_pstd.F90)

set( MAXWELL_SOLVERS_DEPS
  sll_fft
  sll_fields 
  sll_utilities
  sll_linear_solvers
  sll_splines
  sll_fem_spline
  sll_collective
  sll_poisson_solvers
  sll_mapping)

add_library(sll_maxwell_solvers STATIC ${MAXWELL_SOLVERS_SRCS})

target_link_libraries(sll_maxwell_solvers ${MAXWELL_SOLVERS_DEPS})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
