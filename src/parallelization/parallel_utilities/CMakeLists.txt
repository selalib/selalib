add_library(sll_parallel_utilities STATIC sll_m_buffer_loader_utilities.F90)
target_link_libraries(sll_parallel_utilities sll_remap sll_collective)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
