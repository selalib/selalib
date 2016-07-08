ADD_LIBRARY( sll_pic_interface STATIC 
  sll_m_pic_poisson_base.F90
  sll_m_particle_initializer.F90
  sll_m_control_variate.F90
  sll_m_pic_poisson_2d.F90)
TARGET_LINK_LIBRARIES( sll_pic_interface 
  sll_memory
  sll_working_precision
  sll_collective
  sll_poisson_solvers
  sll_maxwell_solvers
  sobol
  prob)


# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
