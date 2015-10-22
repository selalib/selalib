# Sequential library
ADD_LIBRARY( sll_quasi_neutral_solvers STATIC sll_m_qn_solver_2d.F90 )
TARGET_LINK_LIBRARIES( sll_quasi_neutral_solvers sll_fft sll_linear_solvers sll_timer )


# Ctests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)

