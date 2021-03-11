add_library(
  sll_pic_time_integration STATIC
  sll_m_operator_splitting_pic_vp_2d2v.F90
  sll_m_hamiltonian_splitting_pic_vm_1d2v.F90
  sll_m_hamiltonian_splitting_pic_vm_1d2v_boris.F90)

target_link_libraries(
  sll_pic_time_integration
  sll_memory
  sll_assert
  sll_splitting_methods
  sll_particle_groups
  sll_pic_weights
  sll_particle_mesh_coupling
  sll_pic_field_solvers
  sll_poisson_solvers
  sll_maxwell_solvers
  sll_collective)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif()
