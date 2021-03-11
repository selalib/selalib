add_library(sll_descriptors STATIC sll_m_descriptors.F90)

target_link_libraries(sll_descriptors sll_memory)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
