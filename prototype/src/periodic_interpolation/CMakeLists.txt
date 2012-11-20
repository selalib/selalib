IF(FFTPACK_ENABLED)
  add_library(sll_periodic_interp STATIC periodic_interp.F90 bsplvb.f90)
  ADD_DEPENDENCIES(sll_periodic_interp 
    sll_utilities
    sll_fft
    sll_assert
    sll_memory
    sll_working_precision
    sll_splines)
  add_executable(test_periodic_interp unit_test.F90)
  ADD_DEPENDENCIES(test_periodic_interp
    sll_periodic_interp 
    sll_utilities
    sll_fft
    sll_assert
    sll_memory
    sll_working_precision
    sll_splines)
  target_link_libraries(test_periodic_interp sll_periodic_interp
      sll_splines
      sll_fft
      fftpack_module
      sll_utilities
      sll_assertion
      sll_memory
      sll_working_precision
  )

  add_definitions(-D_FFTPACK)
 ENDIF()
