set(SRC_UTILITIES sll_m_utilities.F90 sll_m_fornberg.F90 sll_m_horner.F90)

add_library(sll_utilities STATIC ${SRC_UTILITIES})
target_link_libraries(sll_utilities sll_errors sll_memory sll_assert)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
