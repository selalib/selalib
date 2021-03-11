add_library(sll_random_deviate_generators STATIC sll_m_gaussian.F90
                                                 sll_m_hammersley.F90)
target_link_libraries(sll_random_deviate_generators sll_assert sll_memory)

add_subdirectory(moment_matching)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
