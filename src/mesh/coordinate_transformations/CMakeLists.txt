ADD_LIBRARY( sll_coordinate_transformations STATIC
  sll_m_coordinate_transformation_2d_base.F90
  sll_m_coordinate_transformations_2d.F90
  sll_m_common_coordinate_transformations.F90
  sll_m_singular_mapping_base.F90
  sll_m_singular_mapping_discrete.F90
  sll_m_singular_mapping_analytic.F90
  sll_m_singular_mapping_analytic_target.F90
  sll_m_singular_mapping_analytic_czarny.F90
  sll_m_singular_mapping_advector_base.F90
  sll_m_singular_mapping_advector_constant.F90
  sll_m_singular_mapping_advector_rotating.F90
  sll_m_jacobian_2d_pseudo_cartesian.F90 )

TARGET_LINK_LIBRARIES( sll_coordinate_transformations
  sll_interpolators
  sll_splines
  sll_file_io
  sll_meshes )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
