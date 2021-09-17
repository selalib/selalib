SET ( FEM_SPLINE_SRCS
  sll_m_spline_fem_utilities.F90)
SET( FEM_SPLINE_DEPS
  sll_working_precision
  sll_fft
  sll_assert
  sll_errors
  sll_memory
  sll_constants
  sll_splines
  sll_spline_matrix_helper
  sll_quadrature
  sll_linear_solvers
  sll_distribution_function
  sll_mapping)


set( FEM_SPLINE_SRCS
  ${FEM_SPLINE_SRCS}
  sll_m_spline_fem_utilities_sparse.F90
  sll_m_spline_fem_utilities_3d.F90
  sll_m_spline_fem_utilities_3d_clamped.F90
  sll_m_linear_solver_spline_mass_fft.F90
  sll_m_linear_solver_spline_mass_2d_fft.F90)

set( FEM_SPLINE_DEPS
  ${FEM_SPLINE_DEPS})


ADD_LIBRARY(sll_fem_spline STATIC ${FEM_SPLINE_SRCS})

TARGET_LINK_LIBRARIES(sll_fem_spline ${FEM_SPLINE_DEPS})
