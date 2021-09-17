SET ( SPLINE_MATRIX_HELPER_SRCS
  sll_m_spline_fem_utilities_3d_helper.F90)
SET( SPLINE_MATRIX_HELPER_DEPS
  sll_working_precision
  sll_assert
  sll_memory
  sll_linear_solvers)


set( SPLINE_MATRIX_HELPER_DEPS
  ${SPLINE_MATRIX_HELPER_DEPS})


ADD_LIBRARY(sll_spline_matrix_helper STATIC ${SPLINE_MATRIX_HELPER_SRCS})

TARGET_LINK_LIBRARIES(sll_spline_matrix_helper ${SPLINE_MATRIX_HELPER_DEPS})
