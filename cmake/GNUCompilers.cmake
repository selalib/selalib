# GNU compilers
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions(-Drestrict=__restrict__ -DADD_ -DINLINE_ALL=inline)
  # SET(CMAKE_CXX_FLAGS "-O3 -ftemplate-depth-60 -Drestrict=__restrict__
  # -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math
  # -Wno-deprecated ") SET(CMAKE_CXX_FLAGS "-g -O3 -ftemplate-depth-60
  # -Drestrict=__restrict__ -funroll-all-loops   -finline-limit=1000
  # -Wno-deprecated ")
  set(GNU_FLAGS
      "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated "
  )
  set(CMAKE_CXX_FLAGS "${GNU_FLAGS} -ftemplate-depth=60")
  set(CMAKE_C_FLAGS "${GNU_FLAGS} -std=c99")

  if(HAVE_POSIX_MEMALIGN)
    set(CMAKE_TRY_GNU_CC_FLAGS "-mmmx")
    check_c_compiler_flag(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
    if(GNU_CC_FLAGS)
      set(HAVE_MMX 1)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmmx")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmmx")
    endif(GNU_CC_FLAGS)

    set(CMAKE_TRY_GNU_CC_FLAGS "-msse")
    check_c_compiler_flag(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
    if(GNU_CC_FLAGS)
      set(HAVE_SSE 1)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
    endif(GNU_CC_FLAGS)

    set(CMAKE_TRY_GNU_CXX_FLAGS "-msse2")
    check_c_compiler_flag(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
    if(GNU_CC_FLAGS)
      set(HAVE_SSE2 1)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    endif(GNU_CC_FLAGS)

    set(CMAKE_TRY_GNU_CC_FLAGS "-msse3")
    check_c_compiler_flag(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
    if(GNU_CC_FLAGS)
      set(HAVE_SSE3 1)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3")
    endif(GNU_CC_FLAGS)
  endif(HAVE_POSIX_MEMALIGN)

  # SET(CMAKE_CXX_FLAGS "-O6 -ftemplate-depth-60 -Drestrict=__restrict__
  # -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math
  # -Wno-deprecated -pg") SET(CMAKE_CXX_FLAGS "-g -ftemplate-depth-60
  # -Drestrict=__restrict__ -fstrict-aliasing -Wno-deprecated")

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS}
    # -faltivec -framework Accelerate -bind_at_load")
    set(F77 xlf)
    set(F77FLAGS -O3)
  else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # SET(FORTRAN_LIBS "-lg2c")
    set(F77 g77)
    set(F77FLAGS -funroll-loops -O3)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  # INCLUDE(${CMAKE_ROOT}/Modules/TestCXXAcceptsFlag.cmake) IF(QMC_OMP)
  set(CMAKE_TRY_OPENMP_CXX_FLAGS "-fopenmp")
  check_cxx_accepts_flag(${CMAKE_TRY_OPENMP_CXX_FLAGS} GNU_OPENMP_FLAGS)
  if(GNU_OPENMP_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_TRY_OPENMP_CXX_FLAGS}")
    set(ENABLE_OPENMP 1)
  endif(GNU_OPENMP_FLAGS)
  # ENDIF(QMC_OMP)

  if(QMC_BUILD_STATIC)
    set(CMAKE_CXX_LINK_FLAGS " -static")
  endif(QMC_BUILD_STATIC)

  set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${CMAKE_C_FLAGS}")

endif(CMAKE_COMPILER_IS_GNUCXX)

# IF(APPLE) INCLUDE_DIRECTORIES(/sw/include) ENDIF(APPLE)
