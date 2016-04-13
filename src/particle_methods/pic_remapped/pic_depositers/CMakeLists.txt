ADD_LIBRARY( sll_pic_depositers STATIC
  sll_m_pic_depositers.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_depositers
  #  sll_memory
  sll_assert
  sll_pic_interface
  sll_pic_lbfr
  sll_pic_opt2d_accumulators
  )
