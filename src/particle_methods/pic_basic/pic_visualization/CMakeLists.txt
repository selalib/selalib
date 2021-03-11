add_library(sll_pic_visualization STATIC
            sll_m_particle_visualization_interface.F90)

target_link_libraries(
  sll_pic_visualization
  sll_assert
  sll_errors
  sll_memory
  sll_working_precision
  sll_particle_groups
  sll_meshes
  sll_file_io)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
