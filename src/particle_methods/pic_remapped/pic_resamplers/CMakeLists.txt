ADD_LIBRARY( sll_pic_resamplers STATIC
  sll_m_pic_resamplers.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_resamplers
  #  sll_memory
  sll_assert
  #  sll_meshes
  #  sll_remapped_pic
  )
