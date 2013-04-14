SET(BSL_LIB_DEPS sll_interpolators
                 sll_splines 
                 sll_memory 
	         sll_working_precision
                 sll_utilities 
                 sll_assertion)

ADD_EXECUTABLE(bsl_1d_cubic_uniform_periodic test_1d_cubic_uniform_periodic.F90)
ADD_DEPENDENCIES(bsl_1d_cubic_uniform_periodic ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(bsl_1d_cubic_uniform_periodic ${BSL_LIB_DEPS})

ADD_EXECUTABLE(bsl_1d_cubic_uniform_compact test_1d_cubic_uniform_compact.F90)
ADD_DEPENDENCIES(bsl_1d_cubic_uniform_compact ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(bsl_1d_cubic_uniform_compact ${BSL_LIB_DEPS})

ADD_EXECUTABLE(bsl_1d_cubic_nonuniform_periodic test_1d_cubic_nonuniform_periodic.F90)
ADD_DEPENDENCIES(bsl_1d_cubic_nonuniform_periodic ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(bsl_1d_cubic_nonuniform_periodic ${BSL_LIB_DEPS})

ADD_EXECUTABLE(bsl_1d_cubic_nonuniform_compact test_1d_cubic_nonuniform_compact.F90)
ADD_DEPENDENCIES(bsl_1d_cubic_nonuniform_compact ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(bsl_1d_cubic_nonuniform_compact ${BSL_LIB_DEPS})

ADD_EXECUTABLE(bsl_1d_quintic_uniform_compact test_1d_quintic_uniform_compact.F90)
ADD_DEPENDENCIES(bsl_1d_quintic_uniform_compact ${BSL_LIB_DEPS})
ADD_EXECUTABLE(bsl_1d_quintic_nonuniform_compact test_1d_quintic_nonuniform_compact.F90)
ADD_DEPENDENCIES(bsl_1d_quintic_nonuniform_compact ${BSL_LIB_DEPS})

IF (LAPACK_LIBRARIES)
   TARGET_LINK_LIBRARIES(bsl_1d_quintic_uniform_compact ${BSL_LIB_DEPS} ${LAPACK_LIBRARIES})
   TARGET_LINK_LIBRARIES(bsl_1d_quintic_nonuniform_compact ${BSL_LIB_DEPS} ${LAPACK_LIBRARIES})
ELSE()
   TARGET_LINK_LIBRARIES(bsl_1d_quintic_uniform_compact ${BSL_LIB_DEPS})
   TARGET_LINK_LIBRARIES(bsl_1d_quintic_nonuniform_compact ${BSL_LIB_DEPS})
ENDIF()


IF(NOT STDF95)

   ADD_EXECUTABLE(bsl_2d_omp unit_test_2d_omp.F90)
   ADD_DEPENDENCIES(bsl_2d_omp ${BSL_LIB_DEPS})
   TARGET_LINK_LIBRARIES(bsl_2d_omp ${BSL_LIB_DEPS})

   INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../selalib_level)
   ADD_EXECUTABLE(bsl_landau_1d landau_1d.F90)
   ADD_DEPENDENCIES(bsl_landau_1d selalib)
   TARGET_LINK_LIBRARIES(bsl_landau_1d selalib dfftpack)

   IF(FFTW_FOUND)

      ADD_EXECUTABLE(bsl_landau_2d landau_2d.F90)
      ADD_DEPENDENCIES(bsl_landau_2d ${BSL_LIB_DEPS})
      TARGET_LINK_LIBRARIES(bsl_landau_2d ${BSL_LIB_DEPS} sll_maxwell_solvers ${FFTW_LIBRARIES})

   ENDIF()

ENDIF(NOT STDF95)

#ADD_LIBRARY(sll_bsl STATIC sll_bsl.F90)

