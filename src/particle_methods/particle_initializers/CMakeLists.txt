ADD_LIBRARY( sll_pic_random_initializers STATIC pic_random_initializers.F90 )

TARGET_LINK_LIBRARIES( sll_pic_random_initializers
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes
  sll_pic_interface
  sll_random_deviate_generators
  )
