ADD_LIBRARY( sll_pic_initializers STATIC
  sll_m_pic_initializers.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_initializers
  #  sll_memory
  sll_assert
  sll_pic_interface
  sll_pic_lbfr
  )
