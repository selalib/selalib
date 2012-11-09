SET(BSL_LIB_DEPS sll_interpolators
                 sll_splines 
                 sll_memory 
	         sll_working_precision
                 sll_utilities 
                 sll_assertion)

ADD_EXECUTABLE(test_bsl_1d unit_test_1d.F90)
ADD_DEPENDENCIES(test_bsl_1d ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(test_bsl_1d ${BSL_LIB_DEPS})

ADD_EXECUTABLE(test_bsl_2d unit_test_2d.F90)
ADD_DEPENDENCIES(test_bsl_2d ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(test_bsl_2d ${BSL_LIB_DEPS})

IF(FFTPACK_ENABLED)
   ADD_EXECUTABLE(test_bsl_landau_2d landau_2d.F90)
   ADD_DEPENDENCIES(test_bsl_landau_2d ${BSL_LIB_DEPS})
   TARGET_LINK_LIBRARIES(test_bsl_landau_2d ${BSL_LIB_DEPS} 
                                            sll_poisson_1d_periodic
                                            fftpack_module)
ENDIF()

IF(FFTW_FOUND)
   ADD_EXECUTABLE(test_bsl_landau_4d landau_4d.F90)
   ADD_DEPENDENCIES(test_bsl_landau_4d ${BSL_LIB_DEPS})
   TARGET_LINK_LIBRARIES(test_bsl_landau_4d ${BSL_LIB_DEPS}
                                            sll_maxwell
                                            ${FFTW_LIBRARIES})
ENDIF()

#ADD_LIBRARY(sll_bsl STATIC sll_bsl.F90)

