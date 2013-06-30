SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

ADD_LIBRARY(sll_mgd2 gxch1_2d.F90
                     gscale_2d.F90
                     diagnostics.F90
                     write_xdmf_2d.F90
                     mgd2.F90
                     mgd_polar.F90)

ADD_LIBRARY(sll_mgd3 mgdend_3d.F90 
                     mgdkcyc_3d.F90
                     mgd3.F90
                     gscale_3d.F90
                     gxch1cor_3d.F90
                     gxch1lin_3d.F90
                     gxch1pla_3d.F90
                     mgdbdry_3d.F90 
                     mgdcor_3d.F90 
                     mgderr_3d.F90 
                     mgdpfpde_3d.F90 
                     mgdphpde_3d.F90 
                     mgdppde_3d.F90 
                     mgdrelax_3d.F90
                     mgdrestr_3d.F90
                     mgdrpde_3d.F90 
                     mgdrsetf_3d.F90
                     mgdsetf_3d.F90 
                     mgdrtrsf_3d.F90
                     write_xdmf_3d.F90)

ADD_LIBRARY(sll_multigrid sll_multigrid_base.F90
                          sll_multigrid_2d.F90
                          sll_multigrid_polar.F90)

TARGET_LINK_LIBRARIES(sll_multigrid sll_mgd2 
                                    sll_remap 
                                    sll_utilities
                                    sll_collective 
                                    sll_assertion
                                    sll_memory)

ADD_EXECUTABLE(test_multigrid_2d test_multigrid_2d.F90)

TARGET_LINK_LIBRARIES(test_multigrid_2d sll_multigrid sll_file_io )

IF(HDF5_PARALLEL_ENABLED)

   ADD_EXECUTABLE(test_mgd3 unit_test_3d.F90)

   TARGET_LINK_LIBRARIES(test_mgd3 sll_mgd3 
                                   sll_file_io_parallel
	                             sll_constants
	                             sll_utilities
	                             sll_collective
	                             sll_memory
	                             sll_assertion
	                             sll_working_precision )


   ADD_EXECUTABLE(test_mgd2 unit_test_2d.F90)

   TARGET_LINK_LIBRARIES(test_mgd2 sll_mgd2 
                                   sll_file_io_parallel
	                             sll_collective
	                             sll_utilities
	                             sll_memory
	                             sll_assertion)
ENDIF()
