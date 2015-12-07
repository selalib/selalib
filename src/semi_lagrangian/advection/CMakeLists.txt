IF(MPI_ENABLED)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF()

SET( ADVECTION_SRC_FILES 
       sll_m_lagrange_interpolation.F90 
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
       sll_m_advection_2d_tri_mesh.F90)

ADD_LIBRARY( sll_advection STATIC ${ADVECTION_SRC_FILES})

TARGET_LINK_LIBRARIES( sll_advection
                       sll_meshes
                       sll_characteristics
                       sll_quadrature
                       sll_periodic_interpolation
                       sll_splines
                       sll_fft
                       sll_utilities)

IF(BUILD_TESTING)

  ADD_SUBDIRECTORY(testing)

ENDIF(BUILD_TESTING)
