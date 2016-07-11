ADD_LIBRARY( sll_pic_field_solvers STATIC 
  sll_m_pic_poisson_base.F90
  sll_m_pic_poisson_2d.F90)
TARGET_LINK_LIBRARIES( sll_pic_field_solvers 
  sll_memory
  sll_working_precision
  sll_collective
  sll_poisson_solvers
  sll_maxwell_solvers)


# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
