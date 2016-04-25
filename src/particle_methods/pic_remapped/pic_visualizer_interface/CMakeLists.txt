ADD_LIBRARY( sll_pic_visualizer_interface STATIC
  sll_m_pic_visualizer_interface.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_visualizer_interface
  #  sll_memory
  sll_assert
  sll_pic_interface
  sll_pic_lbfr
  )
