ADD_LIBRARY( sll_initial_density_parameters STATIC
  sll_m_initial_density_parameters.F90
  )

TARGET_LINK_LIBRARIES( sll_initial_density_parameters
  #  sll_memory
  sll_assert
  )

ADD_LIBRARY( sll_particle_sampling_interface STATIC
  sll_m_particle_sampling_interface.F90
  )

TARGET_LINK_LIBRARIES( sll_particle_sampling_interface
  #  sll_memory
  sll_assert
  sll_pic_interface
  sll_initial_density_parameters
  sll_pic_lbfr
  sll_pic_remapped_lbf
  )

