SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
ADD_LIBRARY(selalib_vlasov4d STATIC sll_vlasov4d_spectral.F90
                                    sll_vlasov4d_spectral_charge.F90
                                    sll_vlasov4d_maxwell.F90
                                    sll_vlasov4d_poisson.F90
                                    sll_vlasov4d_polar.F90
                                    sll_vlasov4d_base.F90
                                    init_functions.F90)

TARGET_LINK_LIBRARIES(selalib_vlasov4d sll_assert 
  sll_file_io sll_file_io_parallel sll_remap sll_utilities 
  sll_maxwell_solvers sll_poisson_solvers sll_collective sll_fft sll_mpi)

ADD_EXECUTABLE(vp4d_remapper        vp4d_remapper.F90)
ADD_EXECUTABLE(vm4d_remapper        vm4d_remapper.F90)
ADD_EXECUTABLE(vm4d_spectral        vm4d_spectral.F90)
ADD_EXECUTABLE(vp4d_polar           vp4d_polar.F90)
ADD_EXECUTABLE(vm4d_spectral_charge vm4d_spectral_charge.F90)
ADD_EXECUTABLE(vm2d_spectral_charge vm2d_spectral_charge.F90)
ADD_EXECUTABLE(vp4d_multigrid       vp4d_multigrid.F90)
ADD_EXECUTABLE(vp4d_sequential      vp4d_sequential.F90)

TARGET_LINK_LIBRARIES(vp4d_remapper        selalib_vlasov4d)
TARGET_LINK_LIBRARIES(vm4d_remapper        selalib_vlasov4d)
TARGET_LINK_LIBRARIES(vm4d_spectral        selalib_vlasov4d)
TARGET_LINK_LIBRARIES(vp4d_polar           selalib_vlasov4d)
TARGET_LINK_LIBRARIES(vm4d_spectral_charge selalib_vlasov4d)
TARGET_LINK_LIBRARIES(vm2d_spectral_charge selalib_vlasov4d)
TARGET_LINK_LIBRARIES(vp4d_sequential      sll_fft)
TARGET_LINK_LIBRARIES(vp4d_multigrid       selalib_vlasov4d)
   
FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../tests DESTINATION ${CMAKE_BINARY_DIR})
