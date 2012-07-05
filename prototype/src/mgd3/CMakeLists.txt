set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -fdefault-real-8 -ffree-line-length-none -fall-intrinsics")
add_library(mgd3 STATIC decomp.F gerr.F ginit.F grid1_type.F gscale.F gxch1cor.F
gxch1lin.F gxch1pla.F mgdbdry.F mgdcor.F mgdend.F mgderr.F mgdinit.F mgdkcyc.F
mgdpfpde.F mgdphpde.F mgdppde.F mgdrelax.F mgdrestr.F mgdrpde.F mgdrsetf.F
mgdrtrsf.F mgdsetf.F mgdsolver.F)

add_executable(tmgd3 tmgd3.F)

target_link_libraries(tmgd3 mgd3)
