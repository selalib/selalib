ADD_LIBRARY( sll_pic_sampling STATIC 
  sll_m_particle_initializer.F90)
TARGET_LINK_LIBRARIES( sll_pic_sampling 
  sll_memory
  sll_errors
  sll_working_precision
  sll_particle_groups
  sobol
  prob)


# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
