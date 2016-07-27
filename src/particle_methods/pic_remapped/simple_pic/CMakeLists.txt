ADD_LIBRARY( sll_simple_pic_4d_particle STATIC
  sll_m_simple_pic_4d_particle.F90)

TARGET_LINK_LIBRARIES( sll_simple_pic_4d_particle
  sll_working_precision
  sll_memory)

ADD_LIBRARY( sll_simple_pic STATIC
  sll_m_simple_pic_4d_group.F90
  )

TARGET_LINK_LIBRARIES( sll_simple_pic
  sll_assert
  sll_meshes
  sll_remapped_pic
  sll_particle_groups
  sll_pic_opt2d_accumulators
  sll_remapped_pic_utilities
  sll_simple_pic_4d_particle
  )

# Ctest
IF(BUILD_TESTING)
  add_subdirectory( testing )
ENDIF(BUILD_TESTING)
