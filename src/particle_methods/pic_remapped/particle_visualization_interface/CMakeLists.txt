ADD_LIBRARY( sll_particle_visualization_interface STATIC
  sll_m_particle_visualization_interface.F90
  )

TARGET_LINK_LIBRARIES( sll_particle_visualization_interface
  sll_assert
  sll_pic_interface
  sll_meshes
  sll_pic_lbfr
  sll_pic_remapped_lbf
  sll_file_io
  )
