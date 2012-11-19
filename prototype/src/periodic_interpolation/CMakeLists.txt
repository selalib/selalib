IF(FFTPACK_ENABLED)

  add_definitions(-D_FFTPACK)
  add_library(sll_periodic_interp STATIC periodic_interp.F90 bsplvb.f90)

  ADD_DEPENDENCIES(sll_periodic_interp 
                   sll_splines
                   sll_fft
                   fftpack_module
                   sll_utilities
                   sll_assert
                   sll_memory
                   sll_working_precision)

  add_executable(test_periodic_interp unit_test.F90)

  target_link_libraries(test_periodic_interp sll_periodic_interp
                                             sll_splines
                                             sll_fft
                                             sll_utilities
                                             fftpack_module
                                             sll_assertion
                                             sll_memory
                                             sll_working_precision)

 ENDIF()
