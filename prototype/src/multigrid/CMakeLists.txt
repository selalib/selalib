IF (MPI_MODULE_ENABLED)

  ADD_LIBRARY(mgd2_periodic mgd2_periodic.F)
  ADD_LIBRARY(mgd2_nonperiodic mgd2_nonperiodic.F)

  ADD_EXECUTABLE(test_multigrid_2d_periodic sll_multigrid_2d.F90 
                                            test_multigrid_2d_periodic.F90)

  TARGET_LINK_LIBRARIES(test_multigrid_2d_periodic mgd2_periodic 
                                                   sll_collective 
                                                   sll_utilities
                                                   sll_memory
                                                   sll_working_precision)

  ADD_EXECUTABLE(test_multigrid_2d_dirichlet sll_multigrid_2d.F90 
                                             test_multigrid_2d_dirichlet.F90)

  TARGET_LINK_LIBRARIES(test_multigrid_2d_dirichlet mgd2_nonperiodic 
                                                    sll_collective 
                                                    sll_utilities
                                                    sll_memory
                                                    sll_working_precision)

ENDIF (MPI_MODULE_ENABLED)
