ADD_LIBRARY( sll_simple_pic STATIC
  sll_m_simple_pic_4d_particle.F90
  sll_m_simple_pic_4d_group.F90
  )

TARGET_LINK_LIBRARIES( sll_simple_pic
  sll_memory
  sll_assert
  sll_meshes
  sll_remapped_pic
  # sll_pic_interface
  sll_pic_accumulators
  sll_remapped_pic_utilities
  )

# Ctest
IF(BUILD_TESTING)
  add_subdirectory( testing )
ENDIF(BUILD_TESTING)
