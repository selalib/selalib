set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -ffree-line-length-none -fall-intrinsics")
add_library(mgd3 STATIC 
ginit.F90 
mgdend.F90 
mgdinit.F90 
mgdkcyc.F90
mgdsolver.F90 
mgd3.F90
decomp.F90 
gerr.F90
grid1_type.F90
gscale.F90
gxch1cor.F90
gxch1lin.F90
gxch1pla.F90
mgdbdry.F90 
mgdcor.F90 
mgderr.F90 
mgdpfpde.F90 
mgdphpde.F90 
mgdppde.F90 
mgdrelax.F90
mgdrestr.F90
mgdrpde.F90 
mgdrsetf.F90
mgdsetf.F90 
mgdrtrsf.F90
)

add_executable(tmgd3 tmgd3.F90)

target_link_libraries(tmgd3 mgd3)
