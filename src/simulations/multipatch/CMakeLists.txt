IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)

  ADD_SUBDIRECTORY(distribution_function_multipatch)
  ADD_SUBDIRECTORY(meshes_multipatch)
  ADD_SUBDIRECTORY(coordinate_transformations_multipatch)
  ADD_SUBDIRECTORY(fields_multipatch)
  ADD_SUBDIRECTORY(general_coordinate_elliptic_solvers_multipatch)
  ADD_SUBDIRECTORY(sim_4d_qns_general_multipatch)

ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)
