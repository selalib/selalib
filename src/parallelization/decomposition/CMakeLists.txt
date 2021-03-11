add_library(sll_decomposition STATIC sll_m_decomposition.F90)
target_link_libraries(sll_decomposition sll_collective)

# Ctests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
