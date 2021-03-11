# Set correct compilation flags for FFTPACK, depending on Fortran compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -w -cpp -fdefault-real-8 ")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -w -cpp -fdefault-real-8 ")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -r8 -nowarn")
  set(CMAKE_Fortran_FLAGS_RELEASE "-r8 -nowarn -O3 -xHost -ip -qopenmp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES IBM)
  set(CMAKE_Fortran_COMPILER xlf)
  set(CMAKE_Fortran_FLAGS_DEBUG "-qfixed -qrealsize=8")
  set(CMAKE_Fortran_FLAGS_RELEASE "-qfixed -qrealsize=8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -r8 ")
  set(CMAKE_Fortran_FLAGS_RELEASE "-fast -r8 ")
endif()

# New version of FFTPACK (this will substitute old 'dfftpack')
add_library(fftpack5 STATIC fftpack5.f)

# Every external library should be installed (otherwise 'package' does not work)
install(TARGETS fftpack5 ARCHIVE DESTINATION lib)

# CTest
if(BUILD_TESTING)
  add_executable(test_fftpack5 fftpack5_prb.f)
  target_link_libraries(test_fftpack5 fftpack5)
endif(BUILD_TESTING)
