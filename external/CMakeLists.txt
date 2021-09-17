# Allow compiler-specific intrinsic functions (non-standard)
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fall-intrinsics")
  set(CMAKE_Fortran_FLAGS_RELEASE
      "${CMAKE_Fortran_FLAGS_RELEASE} -fall-intrinsics")
endif()

add_subdirectory(fftpack)
add_subdirectory(pppack)
add_subdirectory(nufft)
add_subdirectory(burkardt)
add_subdirectory(sparsekit)

if(BUILD_MUDPACK)
  add_subdirectory(mudpack)
endif(BUILD_MUDPACK)
