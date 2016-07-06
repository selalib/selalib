ADD_LIBRARY( sll_particle_visualization_interface STATIC
  sll_m_particle_visualization_interface.F90
  )

TARGET_LINK_LIBRARIES( sll_particle_visualization_interface
  sll_assert
  sll_errors
  sll_memory
  sll_working_precision
  sll_pic_interface
  sll_particle_types
  sll_meshes
  sll_file_io
  )
