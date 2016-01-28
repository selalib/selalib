ADD_LIBRARY( sll_ode_solvers STATIC sll_m_ode_solvers.F90 )
TARGET_LINK_LIBRARIES( sll_ode_solvers sll_distribution_function )

# Ctest
IF(BUILD_TESTING)

  ADD_SUBDIRECTORY(testing)

ENDIF(BUILD_TESTING)
