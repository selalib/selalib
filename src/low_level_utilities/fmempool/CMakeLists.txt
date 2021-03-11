add_library(fmempool STATIC fmempool.F90 fassert.F90)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
