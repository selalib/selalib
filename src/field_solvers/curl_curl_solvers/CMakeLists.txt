set(CURL_CURL_SOLVERS_SRCS
  sll_m_preconditioner_jacobi.F90
  sll_m_preconditioner_curl_solver_fft.F90
  sll_m_linear_operator_MG.F90
  sll_m_linear_operator_GTM.F90
  sll_m_linear_operator_curl_3d.F90
  sll_m_linear_operator_MG_cl.F90
  sll_m_linear_operator_GTM_cl.F90
  sll_m_linear_operator_curl_cl_3d.F90)

set( CURL_CURL_SOLVERS_DEPS
  sll_maxwell_solvers
  sll_fft
  sll_utilities
  sll_linear_solvers
  sll_splines
  sll_fem_spline
  sll_mapping)

add_library(sll_curl_curl_solvers STATIC ${CURL_CURL_SOLVERS_SRCS})

target_link_libraries(sll_curl_curl_solvers ${CURL_CURL_SOLVERS_DEPS})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
