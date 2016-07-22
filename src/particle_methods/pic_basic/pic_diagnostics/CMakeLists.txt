ADD_LIBRARY(sll_pic_diagnostics STATIC 
  sll_m_pic_diagnostics.F90)

TARGET_LINK_LIBRARIES( sll_pic_diagnostics
  sll_memory
  sll_particle_groups
  sll_particle_mesh_coupling
  sll_maxwell_solvers
  )


## Tests
#IF(BUILD_TESTING)
#  ADD_SUBDIRECTORY(testing)
#ENDIF()
