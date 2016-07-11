ADD_LIBRARY( sll_initial_densities_temp STATIC
  sll_m_initial_density_parameters.F90
  )

TARGET_LINK_LIBRARIES( sll_initial_densities_temp
  sll_working_precision
  )
