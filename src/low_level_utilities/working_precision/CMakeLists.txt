set(PRECISION_SRC sll_m_working_precision.F90)
add_library(sll_working_precision STATIC ${PRECISION_SRC})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
