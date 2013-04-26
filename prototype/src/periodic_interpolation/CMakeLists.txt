ADD_LIBRARY(sll_periodic_interp STATIC periodic_interp.F90 bsplvb.f90)

ADD_DEPENDENCIES(sll_periodic_interp sll_utilities
                                     sll_fft
                                     sll_assert
                                     sll_memory
                                     sll_working_precision
                                     sll_splines)

ADD_EXECUTABLE(test_periodic_interp unit_test.F90)

TARGET_LINK_LIBRARIES(test_periodic_interp sll_periodic_interp
                                           sll_fft
                                           sll_utilities
                                           sll_splines
                                           dfftpack
                                           sll_assertion
                                           sll_memory
                                           sll_working_precision)
