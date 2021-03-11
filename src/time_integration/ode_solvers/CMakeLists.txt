add_library(sll_ode_solvers STATIC sll_m_ode_solvers.F90)
target_link_libraries(sll_ode_solvers sll_distribution_function)

# Ctest
if(BUILD_TESTING)

  add_subdirectory(testing)

endif(BUILD_TESTING)
