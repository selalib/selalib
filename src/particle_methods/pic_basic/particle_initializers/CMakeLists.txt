ADD_LIBRARY( sll_initial_density_parameters STATIC
  sll_m_initial_density_parameters.F90
  )

TARGET_LINK_LIBRARIES( sll_initial_density_parameters
  sll_working_precision
  )

# ADD_LIBRARY(sll_particle_initializers STATIC sll_m_particle_initializer.F90 ) # old library in this folder but now file has moved

ADD_LIBRARY(sll_particle_sampling_interface STATIC
  sll_m_particle_sampling_interface.F90
  )

TARGET_LINK_LIBRARIES( sll_particle_sampling_interface
  sll_assert
  sll_errors
  sll_memory
  sll_working_precision
  sll_pic_interface                 # for the random particle initializers
  sll_initial_density_parameters
  sll_particle_types
  sobol # not sure this is needed
  prob  # not sure this is needed
  )
