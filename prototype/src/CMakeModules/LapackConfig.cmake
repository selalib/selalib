
IF (APPLE)
   FIND_LIBRARY(LAPACK_LIBRARIES  lapack)
   FIND_LIBRARY(BLAS_LIBRARIES    blas)
ELSE()
   FIND_PACKAGE(LAPACK REQUIRED)
   FIND_PACKAGE(BLAS   REQUIRED)
ENDIF()

SET(LINK_LIBRARIES  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
