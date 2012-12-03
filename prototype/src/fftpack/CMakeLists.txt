set(CMAKE_Fortran_FLAGS "-O3")
add_library(fftpack_module STATIC dfftpack.f fftpack.f90)
SET(CMAKE_Fortran_FLAGS "")
