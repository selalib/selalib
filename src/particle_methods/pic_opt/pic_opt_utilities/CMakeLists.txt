add_library(
  sll_pic_opt_utilities STATIC
  sll_m_pic_utilities.F90 sll_m_charge_to_density.F90 sll_m_particle_sort.F90)

target_link_libraries(
  sll_pic_opt_utilities
  sll_pic_opt_particles
  sll_pic_opt_accumulators
  sll_pic_opt_particles
  sll_assert
  sll_meshes
  sll_memory
  sll_timer
  sll_file_io)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
