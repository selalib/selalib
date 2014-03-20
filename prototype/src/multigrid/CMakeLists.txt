IF (MPI_MODULE_ENABLED)

  ADD_LIBRARY(sll_multigrid_2d 
                  decomp.F 
                  gxch1lin.F 
                  gxch1cor.F 
                  grid1_type.F
		      ginit.F 
                  gscale.F 
                  gerr.F 
		      sll_multigrid_2d.F 
                  diagnostics.f90
             )

  ADD_EXECUTABLE(test_multigrid_2d main.F)

  TARGET_LINK_LIBRARIES(test_multigrid_2d sll_multigrid_2d)

ENDIF (MPI_MODULE_ENABLED)
