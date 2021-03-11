add_library(sll_constants STATIC sll_m_constants.F90)
target_link_libraries(sll_constants sll_working_precision)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
