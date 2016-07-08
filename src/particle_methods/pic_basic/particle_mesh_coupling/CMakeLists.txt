ADD_LIBRARY(sll_particle_mesh_coupling STATIC 
  sll_m_particle_mesh_coupling_base.F90
  sll_m_particle_mesh_coupling_spline_1d.F90
  sll_m_particle_mesh_coupling_spline_2d.F90)

TARGET_LINK_LIBRARIES( sll_particle_mesh_coupling
  sll_memory
  sll_assert
  sll_pic_interface
  sll_particle_groups
  sll_splines
  )


# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
