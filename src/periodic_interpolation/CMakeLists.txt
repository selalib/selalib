ADD_LIBRARY( sll_periodic_interpolation STATIC periodic_interp.F90 bsplvb.f90 )
TARGET_LINK_LIBRARIES( sll_periodic_interpolation sll_fft sll_splines dfftpack )

# Ctest
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
