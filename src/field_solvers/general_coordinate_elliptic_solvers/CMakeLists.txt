add_library(
  sll_general_coordinate_elliptic_solvers STATIC
  sll_m_lobatto_poisson.F90
  sll_m_map_function.F90
  sll_m_lobalap.F90
  sll_m_ellipt_2d_fem_sps_weak_form.F90
  sll_m_ellipt_2d_fem_sps_stencil_new_projector.F90
  sll_m_ellipt_2d_cartesian_gradient.F90)

target_link_libraries(
  sll_general_coordinate_elliptic_solvers
  sll_sparse_matrix
  sll_fields
  sll_interpolators
  sll_coordinate_transformations
  sll_meshes
  sll_fields
  sll_file_io
  sll_quadrature
  sll_utilities
  sll_timer
  sll_assert
  sll_memory)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
