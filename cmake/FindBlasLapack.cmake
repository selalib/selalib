find_library(
  BLAS_LIBRARIES
  NAMES blas
  PATHS ${BLASLAPACK_DIR}
  NO_DEFAULT_PATH)
find_library(
  LAPACK_LIBRARIES
  NAMES lapack
  PATHS ${BLASLAPACK_DIR}
  NO_DEFAULT_PATH)

set(BLAS_LIBRARIES
    ${BLAS_LIBRARIES}
    CACHE STRING "BLAS library")
set(LAPACK_LIBRARIES
    ${LAPACK_LIBRARIES}
    CACHE STRING "LAPACK library")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLASLAPACK DEFAULT_MSG BLAS_LIBRARIES
                                  LAPACK_LIBRARIES)

if(BLASLAPACK_FOUND)

  message(
    STATUS
      "BLAS and LAPACK have been found using BLASLAPACK_DIR=${BLASLAPACK_DIR}")
  message(STATUS "BLAS_LIBRARIES   : ${BLAS_LIBRARIES}")
  message(STATUS "LAPACK_LIBRARIES : ${LAPACK_LIBRARIES}")
  mark_as_advanced(BLAS_LIBRARIES LAPACK_LIBRARIES)

else()

  message(
    FATAL_ERROR
      "BLAS and LAPACK could not be found in BLASLAPACK_DIR=${BLASLAPACK_DIR}")

endif(BLASLAPACK_FOUND)
