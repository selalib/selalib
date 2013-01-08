set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

ADD_LIBRARY(sll_multigrid STATIC mgdend.F90 
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
                                 write_xdmf_3d.F90)

IF(HDF5_PARALLEL_ENABLED)

   add_executable(test_multigrid_3d unit_test_3d.F90)
   target_link_libraries(test_multigrid_3d sll_multigrid 
                                           sll_low_level_file_io_parallel
	                                   sll_constants
	                                   sll_utilities
	                                   sll_collective
	                                   sll_memory
	                                   sll_assertion
	                                   sll_working_precision 
                                           ${HDF5_LIBRARIES})

   add_executable(test_multigrid_2d unit_test_2d.F90)
   target_link_libraries(test_multigrid_2d sll_multigrid 
                                           sll_low_level_file_io_parallel
	                                   sll_constants
	                                   sll_utilities
	                                   sll_collective
	                                   sll_memory
	                                   sll_assertion
	                                   sll_working_precision 
                                           ${HDF5_LIBRARIES})

ENDIF()


