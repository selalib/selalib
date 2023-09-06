add_library(sll_collective STATIC sll_m_collective.F90)
target_link_libraries(sll_collective sll_assert sll_memory 
                      sll_utilities)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
