add_library(fmempool STATIC fmempool.F90 fassert.F90)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
