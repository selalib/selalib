if(HDF5_ENABLED
   AND HDF5_IS_PARALLEL
   AND HDF5_PARALLEL_ENABLED
   AND MPI_ENABLED)

  set(MULTIPATCH_FILES
      sll_m_distribution_function_4d_multipatch.F90
      sll_m_coordinate_transformations_2d_nurbs.F90
      sll_m_coordinate_transformation_multipatch.F90
      sll_m_scalar_field_2d_multipatch.F90
      sll_m_general_coordinate_elliptic_solver_multipatch.F90)

  list(APPEND MULTIPATCH_FILES sll_m_cartesian_meshes_multipatch.F90)

  add_library(sll_multipatch STATIC ${MULTIPATCH_FILES})

  target_link_libraries(
    sll_multipatch
    sll_coordinate_transformations
    sll_fields
    sll_sparse_matrix
    sll_quadrature
    sll_general_coordinate_elliptic_solvers
    sll_remap
    sll_parallel_array_utilities
    sll_timer
    sll_meshes
    sll_assert
    sll_memory)

  if(BUILD_TESTING)
    add_subdirectory(testing)
  endif(BUILD_TESTING)

endif(
  HDF5_ENABLED
  AND HDF5_IS_PARALLEL
  AND HDF5_PARALLEL_ENABLED
  AND MPI_ENABLED)
