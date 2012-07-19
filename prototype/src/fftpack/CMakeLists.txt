set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -ffree-line-length-none -fall-intrinsics")
add_library(fftpack_module STATIC dfftpack.f fftpack.f90)
