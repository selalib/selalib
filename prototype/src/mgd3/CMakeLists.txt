set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -ffree-line-length-none -fall-intrinsics")
add_library(mgd3 STATIC decomp.F gerr.F ginit.F90 grid1_type.F gscale.F gxch1cor.F
gxch1lin.F gxch1pla.F mgdbdry.F mgdcor.F mgdend.F90 mgderr.F mgdinit.F90 mgdkcyc.F90
mgdpfpde.F mgdphpde.F mgdppde.F mgdrelax.F mgdrestr.F mgdrpde.F mgdrsetf.F
mgdrtrsf.F mgdsetf.F mgdsolver.F90 mgd3.F90)

add_executable(tmgd3 tmgd3.F90)

target_link_libraries(tmgd3 mgd3)
