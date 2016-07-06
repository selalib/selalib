ADD_LIBRARY(sll_particle_groups STATIC
  sll_m_particle_group_2d2v.F90 
  sll_m_particle_group_1d2v.F90)

TARGET_LINK_LIBRARIES( sll_particle_groups
  sll_memory
  sll_assert
  sll_pic_interface
  )

IF (BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()