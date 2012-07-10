set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -ffree-line-length-none -fall-intrinsics")
add_library(mgd3 STATIC 
ginit.F90 
mgdend.F90 
mgdinit.F90 
mgdkcyc.F90
mgdsolver.F90 
mgd3.F90
decomp.F 
gerr.F 
grid1_type.F 
gscale.F 
gxch1cor.F
gxch1lin.F 
gxch1pla.F 
mgdbdry.F 
mgdcor.F 
mgderr.F 
mgdpfpde.F 
mgdphpde.F 
mgdppde.F 
mgdrelax.F 
mgdrestr.F 
mgdrpde.F 
mgdrsetf.F
mgdrtrsf.F 
mgdsetf.F 
)

add_executable(tmgd3 tmgd3.F90)

target_link_libraries(tmgd3 mgd3)
