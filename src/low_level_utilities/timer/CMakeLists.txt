add_library(sll_timer STATIC sll_m_timer.F90)
target_link_libraries(sll_timer sll_working_precision)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
