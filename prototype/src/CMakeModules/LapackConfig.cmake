
MESSAGE(STATUS "LAPACK:CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")

IF(CMAKE_Fortran_COMPILER MATCHES "ifort")

   SET(BLA_VENDOR "Intel")

   IF($ENV{MKLROOT} MATCHES "composer")
      MESSAGE(STATUS "MKLROOT:$ENV{MKLROOT}")
      INCLUDE_DIRECTORIES($ENV{MKLROOT}/include)
      IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
         SET(LAPACK_LIBRARIES -L$ENV{MKLROOT}/lib/intel64 -mkl=sequential)
      ELSE()
         SET(LAPACK_LIBRARIES -L$ENV{MKLROOT}/lib/ia32 -mkl=sequential)
      ENDIF()
      SET(LAPACK_FOUND TRUE)
      SET(BLAS_FOUND TRUE)
      SET(BLAS_LIBRARIES  " ")

   ELSEIF($ENV{HOSTNAME} MATCHES "hydra") 
                                      
      INCLUDE_DIRECTORIES($ENV{MKLROOT}/include/intel64/lp64 $ENV{MKLROOT}/include)
      SET(BLAS_LIBRARIES  "-L$ENV{MKLROOT}/lib/intel64 -lmkl_blas95_lp64")
      SET(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_lapack95_lp64 -lmkl_rt -lpthread -lm")
      SET(LAPACK_FOUND TRUE)
      SET(BLAS_FOUND TRUE)

   ELSEIF($ENV{HOSTNAME} MATCHES "hpc-f0*")

      SET(MKLPATH  "/opt/intel/Compiler/11.1/072/mkl/lib/em64t")
      SET(LAPACK_LIBRARIES -L${MKLPATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread)
      SET(LAPACK_FOUND TRUE)
      SET(BLAS_FOUND TRUE)
      SET(BLAS_LIBRARIES  " ")

   ENDIF()

ELSE()

   FIND_PACKAGE(BLAS   QUIET)
   FIND_PACKAGE(LAPACK QUIET)

ENDIF()

IF(LAPACK_FOUND AND BLAS_FOUND)

  MESSAGE(STATUS "LAPACK and BLAS libraries are ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}")

ELSE(LAPACK_FOUND AND BLAS_FOUND)

  MESSAGE(STATUS "Failed to link LAPACK, BLAS, ATLAS libraries with environments")
  MESSAGE(STATUS "Going to search LAPACK standard paths.")
  FIND_LIBRARY(LAPACK_LIBRARIES lapack)
  FIND_LIBRARY(BLAS_LIBRARIES blas)

  IF(LAPACK_LIBRARIES AND BLAS_LIBRARIES)

    MESSAGE(STATUS "LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
    MESSAGE(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
    SET(LAPACK_FOUND TRUE)
    SET(BLAS_FOUND TRUE)

  ELSE()

     MESSAGE(SEND_ERROR "LAPACK NOT FOUND")
     SET(LAPACK_LIBRARIES " ")
     SET(BLAS_LIBRARIES " ")
 
  ENDIF(LAPACK_LIBRARIES AND BLAS_LIBRARIES)

ENDIF(LAPACK_FOUND AND BLAS_FOUND)

MARK_AS_ADVANCED( LAPACK_LIBRARIES
                  BLAS_LIBRARIES    )

SET(LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
