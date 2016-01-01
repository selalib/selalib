IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)


   ADD_LIBRARY(sll_multipatch
     STATIC sll_m_distribution_function_4d_multipatch.F90
            sll_m_cartesian_meshes_multipatch.F90
            sll_m_coordinate_transformation_multipatch.F90
            sll_m_scalar_field_2d_multipatch.F90
            #connectivity_qn.f90
            sll_m_general_coordinate_elliptic_solver_multipatch.F90)

   TARGET_LINK_LIBRARIES(sll_multipatch
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


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)


ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED)
