ADD_LIBRARY(sll_particle_groups STATIC
  sll_m_particle_group_base.F90
  sll_m_particle_group_2d2v.F90 
  sll_m_particle_group_1d2v.F90
  sll_m_particle_lbf_utilities.F90
  sll_m_particle_group_2d2v_lbf.F90
  )

TARGET_LINK_LIBRARIES( sll_particle_groups
  sll_assert
  sll_errors
  sll_memory
  sll_working_precision
  sll_constants
  # sll_pic_interface
  sll_meshes
  sll_sparse_grid
  )

IF (BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()
