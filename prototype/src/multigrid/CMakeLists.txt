set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
#set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -ffree-line-length-none -fall-intrinsics")
include_directories(${HDF5_INCLUDE_DIRS})
add_library(multigrid STATIC 
mgdend.F90 
mgdkcyc.F90
mgd3.F90
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
write_xdmf_3d.F90
)
add_executable(test_multigrid unit_test.F90)
target_link_libraries(test_multigrid multigrid sll_utilities ${HDF5_LIBRARIES})
