add_library( sll_pic_utilities STATIC
             sll_m_particle_method_descriptors.F90 )

target_link_libraries( sll_pic_utilities
  sll_working_precision )