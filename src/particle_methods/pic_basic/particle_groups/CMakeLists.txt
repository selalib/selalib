SET( PARTICLE_GROUP_SRCS
  sll_m_particle_group_base.F90
  sll_m_particle_group_3d1v.F90
  sll_m_particle_group_2d2v.F90 
  sll_m_particle_group_1d2v.F90
  sll_m_particle_group_1d1v.F90)

SET( PARTICLE_GROUP_DEPS
   sll_assert
   sll_errors
   sll_memory
   sll_working_precision
   sll_constants
   sll_meshes
   )

if (NOT (FFT_LIB STREQUAL "SLLFFT"))
   if(PGI)
     LIST( APPEND PARTICLE_GROUP_SRCS sll_m_particle_group_2d2v_lbf_pgi.F90)
   else(PGI)
     LIST( APPEND PARTICLE_GROUP_SRCS sll_m_particle_group_2d2v_lbf.F90)
   endif(PGI)
   LIST( APPEND PARTICLE_GROUP_DEPS sll_distribution_function sll_sparse_grid)
endif()   

ADD_LIBRARY(sll_particle_groups STATIC
  ${PARTICLE_GROUP_SRCS}
  )

TARGET_LINK_LIBRARIES( sll_particle_groups
  ${PARTICLE_GROUP_DEPS}
  )

IF (BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()
