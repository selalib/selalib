# <<<<BSL-LT-PIC_Library>>>> [[file:.]]

ADD_LIBRARY( sll_bsl_lt_pic_4d_particle STATIC sll_m_bsl_lt_pic_4d_particle.F90)
TARGET_LINK_LIBRARIES( sll_bsl_lt_pic_4d_particle sll_working_precision sll_memory)

ADD_LIBRARY( sll_bsl_lt_pic STATIC
  sll_m_bsl_lt_pic_4d_group.F90
  sll_m_bsl_lt_pic_4d_utilities.F90
  )

TARGET_LINK_LIBRARIES( sll_bsl_lt_pic
  sll_working_precision
  sll_memory
  sll_assert
  sll_meshes
  sll_pic_interface
  sll_bsl_lt_pic_4d_particle
  sll_remapped_pic
  sll_remapped_pic_utilities
  # sll_pic_accumulators
  sll_sparse_grid
  sll_pic_opt2d_accumulators
  sll_timer # needed?
  sll_file_io # needed?
  sobol
  )

# Ctest
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )
