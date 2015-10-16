ADD_LIBRARY( sll_general_coordinate_elliptic_solvers STATIC
             sll_m_lobatto_poisson.F90 
             sll_m_map_function.F90
             sll_m_lobalap.F90
             sll_m_gces_dirichlet.F90
             sll_m_gces_dirichlet_homogeneous.F90
             sll_m_gces_full_periodic.F90
             sll_m_general_coordinate_elliptic_solver.F90 )

TARGET_LINK_LIBRARIES( sll_general_coordinate_elliptic_solvers
                       sll_sparse_matrix 
                       sll_fields
                       sll_interpolators
                       sll_coordinate_transformations
                       sll_meshes
                       sll_boundary_condition_descriptors
                       sll_dg_fields
                       sll_file_io
                       sll_integration
                       sll_constants
                       sll_utilities
                       sll_timer
                       sll_assert
                       sll_memory
                       sll_working_precision)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
