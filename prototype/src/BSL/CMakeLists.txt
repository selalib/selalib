SET(BSL_LIB_DEPS sll_interpolators
                 sll_splines 
                 sll_numerical_utilities	
                 sll_memory 
	         sll_working_precision
                 sll_utilities 
                 sll_assertion)


ADD_EXECUTABLE(bsl_1d_periodic unit_test_1d_periodic.F90)
ADD_DEPENDENCIES(bsl_1d_periodic ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(bsl_1d_periodic ${BSL_LIB_DEPS})

ADD_EXECUTABLE(bsl_1d_compact unit_test_1d_compact.F90)
ADD_DEPENDENCIES(bsl_1d_compact ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(bsl_1d_compact ${BSL_LIB_DEPS})

IF(NOT STDF95)

   ADD_EXECUTABLE(bsl_2d_omp unit_test_2d_omp.F90)
   ADD_DEPENDENCIES(bsl_2d_omp ${BSL_LIB_DEPS})
   TARGET_LINK_LIBRARIES(bsl_2d_omp ${BSL_LIB_DEPS})

   IF(FFTPACK_ENABLED)

      INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../poisson_solvers)
      ADD_EXECUTABLE(bsl_landau_1d landau_1d.F90)
      ADD_DEPENDENCIES(bsl_landau_1d ${BSL_LIB_DEPS})
      TARGET_LINK_LIBRARIES(bsl_landau_1d ${BSL_LIB_DEPS} sll_poisson_1d_periodic fftpack_module)

   ENDIF()

   IF(FFTW_FOUND)

      ADD_EXECUTABLE(bsl_landau_2d landau_2d.F90)
      ADD_DEPENDENCIES(bsl_landau_2d ${BSL_LIB_DEPS})
      TARGET_LINK_LIBRARIES(bsl_landau_2d ${BSL_LIB_DEPS} sll_maxwell ${FFTW_LIBRARIES})

   ENDIF()

ENDIF(NOT STDF95)

#ADD_LIBRARY(sll_bsl STATIC sll_bsl.F90)

