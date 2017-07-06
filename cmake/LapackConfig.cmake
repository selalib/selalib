GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)

IF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   IF(DEFINED ENV{MKLROOT})
      IF(USE_MKL_WAS_ON AND NOT USE_MKL)
         MESSAGE(WARNING "when MKLROOT is defined MKL library must be used. Now setting 'USE_MKL=ON'...")
      ENDIF()
      SET(USE_MKL ON CACHE BOOL "Using Intel Math Kernel Library" FORCE)
      SET(USE_MKL_WAS_ON true CACHE INTERNAL "Previous value of USE_MKL flag")
   ELSE()
      MESSAGE(STATUS "Environment variable is not set, please load mkl vars")
   ENDIF()
ENDIF()


IF(USE_MKL)

   SET(BLA_VENDOR "Intel")

   STRING(REGEX REPLACE "^([^:]*):" " " MKLROOT $ENV{MKLROOT})
   MESSAGE(STATUS "MKLROOT:${MKLROOT}")
   INCLUDE_DIRECTORIES(${MKLROOT}/include/intel64/lp64 ${MKLROOT}/include)
   IF(APPLE)
     SET(LAPACK_LIBRARIES "-mkl")
   ELSE()
     # --- NOTE: Recent versions (?>=11.0) Linux ifort support "-mkl" as well
     IF(Fortran_COMPILER_NAME MATCHES "ifort")
        IF(OPENMP_ENABLED)
            SET(LAPACK_LIBRARIES "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl")
        ELSE()
            SET(LAPACK_LIBRARIES "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl")
        ENDIF()
     ELSEIF(Fortran_COMPILER_NAME MATCHES "gfortran")
        IF(OPENMP_ENABLED)
            SET(LAPACK_LIBRARIES "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl")
        ELSE()
            SET(LAPACK_LIBRARIES " -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl")
        ENDIF()
     ELSE()
       MESSAGE(FATAL_ERROR "Don't know how to link MKL using the present compiler.  Intel ifort and GNU gfortran are supported.")
     ENDIF()
   ENDIF(APPLE)
   SET(BLAS_LIBRARIES ${LAPACK_LIBRARIES})
   SET(BLAS_FOUND TRUE)
   SET(LAPACK_FOUND TRUE)

ELSE()
   
   # --- give the fast OpenBLAS library a try
   FIND_LIBRARY(OPENBLAS_LIBRARIES
                openblas
                HINTS
                   ENV OPENBLAS_ROOT  # prioritize custom installation location
                   CMAKE_SYSTEM_LIBRARY_PATH  # also search the default location
                PATH_SUFFIXES
                   lib64 lib
                DOC "OpenBLAS, the free high performance BLAS and LAPACK implementation")
                
   IF (${OPENBLAS_LIBRARIES} MATCHES "OPENBLAS_LIBRARIES-NOTFOUND")
      # --- fall-back to slow reference implementations
      IF(APPLE)
         SET(BLA_VENDOR "Apple")
      ENDIF(APPLE)
      FIND_PACKAGE(BLAS)
      FIND_PACKAGE(LAPACK)
   ELSE()
      SET(BLA_VENDOR "OpenBLAS")
      SET(LAPACK_LIBRARIES ${OPENBLAS_LIBRARIES})
      SET(BLAS_LIBRARIES ${LAPACK_LIBRARIES})
      SET(BLAS_FOUND TRUE)
      SET(LAPACK_FOUND TRUE)
   ENDIF()

ENDIF()

IF(NOT LAPACK_FOUND AND NOT BLAS_FOUND)
  # --- nothing to do here
  
ELSE()

  MESSAGE(STATUS "Failed to link LAPACK, BLAS, OpenBLAS or ATLAS libraries based on the environment information")
  MESSAGE(STATUS "Going to search further in standard paths")
  FIND_LIBRARY(BLAS_LIBRARIES NAMES blas HINTS /opt/local /usr/local PATH_SUFFIXES lib)
  FIND_LIBRARY(LAPACK_LIBRARIES NAMES lapack HINTS /opt/local /usr/local PATH_SUFFIXES lib)

ENDIF(NOT LAPACK_FOUND AND NOT BLAS_FOUND)


IF(LAPACK_LIBRARIES AND BLAS_LIBRARIES)

  MESSAGE(STATUS "BLAS_LIBRARIES:${BLAS_LIBRARIES}")
  MESSAGE(STATUS "LAPACK_LIBRARIES:${LAPACK_LIBRARIES}")

ELSE()

   MESSAGE(SEND_ERROR "LAPACK NOT FOUND")

ENDIF(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
