ADD_LIBRARY( sll_pic_visualizers STATIC
  sll_m_pic_visualizers.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_visualizers
  #  sll_memory
  sll_assert
  sll_pic_interface
  sll_pic_lbfr
  )
