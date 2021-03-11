add_library(sll_splines_smooth STATIC sll_m_polar_bsplines_2d.F90)

target_link_libraries(
  sll_splines_smooth
  sll_splines
  sll_coordinate_transformations
  sll_working_precision
  sll_assert
  sll_errors
  sll_constants
  sll_utilities
  sll_linear_solvers
  sll_meshes
  sll_timer
  sll_quadrature
  pppack
  ${LAPACK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
