ADD_LIBRARY( sll_pic_random_initializer STATIC pic_random_initializer.F90 )

TARGET_LINK_LIBRARIES( sll_pic_random_initializer
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes
  sll_pic_interface )
