## (MCP) not the best name, but there is already a library called name sll_pic_initializers...
ADD_LIBRARY( sll_pic_interface_initializers STATIC pic_random_initializers.F90 )

        # bsl_lt_pic_4d_initializers.F90)

TARGET_LINK_LIBRARIES( sll_pic_interface_initializers
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes
  sll_pic_interface
  sll_bsl_lt_pic_particles
  sll_random_deviate_generators
  sll_bsl_lt_pic_utilities
  )
