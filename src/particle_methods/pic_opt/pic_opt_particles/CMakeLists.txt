add_library(
  sll_pic_opt_particles STATIC
  sll_m_particle_representations.F90 sll_m_particle_group_2d.F90
  sll_m_particle_group_4d.F90 sll_m_particle_group_6d.F90)

target_link_libraries(sll_pic_opt_particles sll_working_precision sll_memory
                      sll_assert sll_meshes)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
