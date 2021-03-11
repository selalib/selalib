add_library(sll_pic_field_solvers STATIC sll_m_pic_poisson_base.F90
                                         sll_m_pic_poisson_2d.F90)
target_link_libraries(
  sll_pic_field_solvers
  sll_memory
  sll_working_precision
  sll_collective
  sll_particle_mesh_coupling
  sll_poisson_solvers
  sll_maxwell_solvers)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
