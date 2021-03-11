add_library(sll_parallel_array_utilities STATIC
            sll_m_parallel_array_initializer.F90)

target_link_libraries(
  sll_parallel_array_utilities
  sll_remap
  sll_coordinate_transformations
  sll_file_io_parallel
  sll_file_io
  sll_utilities
  sll_meshes
  sll_initial_profiles)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
