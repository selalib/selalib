get_filename_component(Fortran_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  if(DEFINED ENV{MKLROOT})
    if(USE_MKL_WAS_ON AND NOT USE_MKL)
      message(
        WARNING
          "when MKLROOT is defined MKL library must be used. Now setting 'USE_MKL=ON'..."
      )
    endif()
    set(USE_MKL
        ON
        CACHE BOOL "Using Intel Math Kernel Library" FORCE)
    set(USE_MKL_WAS_ON
        true
        CACHE INTERNAL "Previous value of USE_MKL flag")
  else()
    message(STATUS "Environment variable is not set, please load mkl vars")
  endif()
endif()

if(USE_MKL)

  set(BLA_VENDOR "Intel")

  string(REGEX REPLACE "^([^:]*):" " " MKLROOT $ENV{MKLROOT})
  message(STATUS "MKLROOT:${MKLROOT}")
  include_directories(${MKLROOT}/include/intel64/lp64 ${MKLROOT}/include)
  if(APPLE)
    set(LAPACK_LIBRARIES "-mkl")
  else()
    # --- NOTE: Recent versions (?>=11.0) Linux ifort support "-mkl" as well
    if(Fortran_COMPILER_NAME MATCHES "ifort")
      if(OPENMP_ENABLED)
        set(LAPACK_LIBRARIES
            "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl"
        )
      else()
        set(LAPACK_LIBRARIES
            "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl"
        )
      endif()
    elseif(Fortran_COMPILER_NAME MATCHES "gfortran")
      if(OPENMP_ENABLED)
        set(LAPACK_LIBRARIES
            "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl"
        )
      else()
        set(LAPACK_LIBRARIES
            " -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl"
        )
      endif()
    else()
      message(
        FATAL_ERROR
          "Don't know how to link MKL using the present compiler.  Intel ifort and GNU gfortran are supported."
      )
    endif()
  endif(APPLE)
  set(BLAS_LIBRARIES ${LAPACK_LIBRARIES})
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)

else()

  # --- give the fast OpenBLAS library a first try with cmake FindLapack
  if(CMAKE_VERSION VERSION_GREATER 3.6.2)
    set(BLA_VENDOR OpenBLAS)
    find_package(BLAS)
    find_package(LAPACK)
  endif()

  if(NOT LAPACK_FOUND)

    # --- give the fast OpenBLAS library a second try
    find_library(
      OPENBLAS_LIBRARIES openblas
      HINTS ENV OPENBLAS_ROOT # prioritize custom installation location
            CMAKE_SYSTEM_LIBRARY_PATH # also search the default location
            /usr/local/opt/openblas # Library location on mac with homebrew
      PATH_SUFFIXES lib64 lib
      DOC "OpenBLAS, the free high performance BLAS and LAPACK implementation")

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(OPENBLAS DEFAULT_MSG OPENBLAS_LIBRARIES)
    include(CheckFortranFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES ${OPENBLAS_LIBRARIES})
    check_fortran_function_exists(dpbtrs OPENBLAS_HAS_LAPACK)

    if(OPENBLAS_FOUND AND OPENBLAS_HAS_LAPACK)
      message(STATUS "OpenBLAS found and has lapack")
      set(BLA_VENDOR OpenBLAS)
      set(LAPACK_LIBRARIES ${OPENBLAS_LIBRARIES})
      set(BLAS_LIBRARIES ${LAPACK_LIBRARIES})
      set(BLAS_FOUND TRUE)
      set(LAPACK_FOUND TRUE)
    else()
      # --- fall-back to slow reference implementations
      if(APPLE)
        set(BLA_VENDOR Apple)
      else()
        unset(BLA_VENDOR)
      endif(APPLE)
      find_package(BLAS)
      find_package(LAPACK)
    endif()

  endif()

endif()

if(NOT LAPACK_FOUND AND NOT BLAS_FOUND)
  # --- nothing to do here

  message(
    STATUS
      "Failed to link LAPACK, BLAS, OpenBLAS or ATLAS libraries based on the environment information"
  )
  message(STATUS "Going to search further in standard paths")
  find_library(
    BLAS_LIBRARIES
    NAMES blas
    HINTS /opt/local /usr/local
    PATH_SUFFIXES lib)
  find_library(
    LAPACK_LIBRARIES
    NAMES lapack
    HINTS /opt/local /usr/local
    PATH_SUFFIXES lib)

endif(NOT LAPACK_FOUND AND NOT BLAS_FOUND)

if(LAPACK_LIBRARIES AND BLAS_LIBRARIES)

  message(STATUS "BLA_VENDOR:${BLA_VENDOR}")
  message(STATUS "BLAS_LIBRARIES:${BLAS_LIBRARIES}")
  message(STATUS "LAPACK_LIBRARIES:${LAPACK_LIBRARIES}")

else()

  message(SEND_ERROR "LAPACK NOT FOUND")

endif(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
