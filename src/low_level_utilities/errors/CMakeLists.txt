add_library(sll_errors STATIC sll_m_errors.F90)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
