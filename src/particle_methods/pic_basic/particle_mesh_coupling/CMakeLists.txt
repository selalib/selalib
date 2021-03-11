add_library(
  sll_particle_mesh_coupling STATIC
  sll_m_particle_mesh_coupling_base.F90
  sll_m_particle_mesh_coupling_spline_1d.F90
  sll_m_particle_mesh_coupling_spline_2d.F90)

target_link_libraries(sll_particle_mesh_coupling sll_memory sll_assert
                      sll_particle_groups sll_splines sll_quadrature)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
