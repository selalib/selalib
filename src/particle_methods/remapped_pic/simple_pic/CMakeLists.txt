ADD_LIBRARY( sll_simple_pic STATIC
  simple_pic_4d_particle.F90
  simple_pic_4d_group.F90
  )

TARGET_LINK_LIBRARIES( sll_simple_pic
  sll_memory
  sll_working_precision
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
