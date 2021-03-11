add_library(sll_periodic_interpolation STATIC sll_m_periodic_interp.F90)
target_link_libraries(sll_periodic_interpolation sll_fft sll_splines pppack
                      dfftpack)

# Ctest
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
