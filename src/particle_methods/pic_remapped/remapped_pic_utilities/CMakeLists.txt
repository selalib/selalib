ADD_LIBRARY( sll_remapped_pic_utilities STATIC
  sll_m_remapped_pic_utilities.F90
  sll_m_pic_random_initializers.F90
  )

TARGET_LINK_LIBRARIES( sll_remapped_pic_utilities
  sll_memory
  sll_assert
  sll_meshes
  sll_remapped_pic
  sll_random_deviate_generators
  )
