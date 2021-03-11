add_library(sll_remap STATIC sll_m_remapper.F90)

target_link_libraries(sll_remap sll_collective sll_utilities)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
