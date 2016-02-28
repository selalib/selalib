ADD_LIBRARY( sll_pic_interface STATIC 
  sll_m_particle_group_base.F90 
  sll_m_kernel_smoother_base.F90)
TARGET_LINK_LIBRARIES( sll_pic_interface sll_memory sll_working_precision)
