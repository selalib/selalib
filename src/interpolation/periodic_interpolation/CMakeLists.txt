ADD_LIBRARY( sll_periodic_interpolation STATIC sll_m_periodic_interp.F90 )
TARGET_LINK_LIBRARIES( sll_periodic_interpolation sll_fft sll_splines pppack dfftpack )

# Ctest
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
