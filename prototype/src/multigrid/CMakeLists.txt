IF (MPI_MODULE_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

  ADD_LIBRARY(mgd2_periodic mgd2_periodic.F)
  ADD_LIBRARY(mgd2_nonperiodic mgd2_nonperiodic.F)
  ADD_LIBRARY(sll_multigrid sll_multigrid_2d.F90 gp_plot2d.F90)

  ADD_EXECUTABLE(test_multigrid_2d_periodic test_multigrid_2d_periodic.F90)

  TARGET_LINK_LIBRARIES(test_multigrid_2d_periodic sll_multigrid
                                                   mgd2_periodic
                                                   sll_file_io
                                                   sll_collective
                                                   )

  ADD_EXECUTABLE(test_multigrid_2d_dirichlet test_multigrid_2d_dirichlet.F90)

  TARGET_LINK_LIBRARIES(test_multigrid_2d_dirichlet sll_multigrid
                                                    mgd2_nonperiodic 
                                                    sll_file_io
                                                    sll_collective
                                                    )

ENDIF (MPI_MODULE_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)
