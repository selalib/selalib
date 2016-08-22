ADD_LIBRARY( sll_pic_visualization STATIC
  sll_m_particle_visualization_interface.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_visualization
  sll_assert
  sll_errors
  sll_memory
  sll_working_precision
  sll_particle_groups
  sll_meshes
  sll_file_io
  )

# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
