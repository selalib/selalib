IF(HDF5_ENABLED)
#hdf5 is parallel and inside aligned_translation_2d.F90
#hdf5 is used
#so the compiler should be parallel
#should be changed
 
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
                       sll_boundary_condition_descriptors
                       sll_utilities)

IF(BUILD_TESTING)

  ADD_SUBDIRECTORY(testing)



  IF(MPI_ENABLED)

     #ADD_EXECUTABLE(
     #  test_aligned_translation_2d_spaghetti
     #  aligned_translation_2d_spaghetti.F90)

     #TARGET_LINK_LIBRARIES(
     #  test_aligned_translation_2d_spaghetti
     #  sll_file_io_parallel
     #  sll_advection
       #sll_fcisl
     #  sll_utilities)

     ADD_EXECUTABLE(
       test_aligned_translation_2d
       aligned_translation_2d.F90)

     TARGET_LINK_LIBRARIES(
       test_aligned_translation_2d
       sll_file_io
       sll_advection
       sll_utilities
       sll_fcisl_toroidal)

  ENDIF(MPI_ENABLED)

  ADD_EXECUTABLE( test_rotation_2d rotation_2d.F90)
  TARGET_LINK_LIBRARIES( test_rotation_2d sll_advection)


ENDIF(BUILD_TESTING)

ENDIF(HDF5_ENABLED)
