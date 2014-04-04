IF (MPI_MODULE_ENABLED)

  ADD_LIBRARY(sll_multigrid_2d 
		      mgd2d_routines.F 
		      sll_multigrid_2d.F90 
                  diagnostics.f90
             )

  ADD_EXECUTABLE(test_multigrid_2d test_multigrid_2d.F90)

  TARGET_LINK_LIBRARIES(test_multigrid_2d sll_multigrid_2d 
                                          sll_utilities
                                          sll_memory
                                          sll_working_precision)

ENDIF (MPI_MODULE_ENABLED)
