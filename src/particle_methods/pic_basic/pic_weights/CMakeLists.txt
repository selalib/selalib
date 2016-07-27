ADD_LIBRARY( sll_pic_weights STATIC 
  sll_m_control_variate.F90)
TARGET_LINK_LIBRARIES( sll_pic_weights
  sll_memory
  sll_working_precision
  sll_particle_groups
  sll_distribution_function)


# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
