add_library(sll_moment_matching STATIC sll_m_moment_matching.F90)
target_link_libraries(sll_moment_matching sll_assert sll_memory)
# Test
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
