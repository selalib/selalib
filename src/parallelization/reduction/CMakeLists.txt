add_library(sll_reduction STATIC sll_m_reduction.F90)
target_link_libraries(sll_reduction sll_meshes)

# Ctest
if(BUILD_TESTING)

  add_subdirectory(testing)

endif(BUILD_TESTING)
