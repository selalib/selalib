add_library(sll_assert STATIC sll_m_assert.F90)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
