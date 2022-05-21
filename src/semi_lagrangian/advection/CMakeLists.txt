set(ADVECTION_SRC_FILES
    sll_m_advection_1d_periodic.F90
    sll_m_advection_1d_non_uniform_cubic_splines.F90
    sll_m_advection_1d_BSL.F90
    sll_m_advection_1d_CSL.F90
    sll_m_advection_1d_CSL_periodic.F90
    sll_m_advection_1d_PSM.F90
    sll_m_advection_1d_base.F90
    sll_m_advection_2d_bsl.F90
    sll_m_advection_2d_integer_oblic.F90
    sll_m_advection_2d_oblic.F90
    sll_m_advection_2d_tensor_product.F90
    sll_m_advection_2d_CSL.F90
    sll_m_advection_2d_base.F90
    sll_m_advection_1d_spectral.F90
    sll_m_advection_1d_ampere.F90
    sll_m_advection_2d_tri_mesh.F90
    sll_m_advection_6d_lagrange_dd_slim.F90
    sll_m_advection_6d_spline_dd_slim.F90)

add_library(sll_advection STATIC ${ADVECTION_SRC_FILES})

target_link_libraries(
  sll_advection
  sll_meshes
  sll_characteristics
  sll_quadrature
  sll_periodic_interpolation
  sll_splines
  sll_fft
  sll_utilities)

if(BUILD_TESTING)

  add_subdirectory(testing)

endif(BUILD_TESTING)
