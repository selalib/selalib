IF( NOT STDF95)
IF(HDF5_ENABLED)
#hdf5 is parallel and inside aligned_translation_2d.F90
#hdf5 is used
#so the compiler should be parallel
#should be changed
 
  IF(MPI_ENABLED)
    SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  ENDIF()

  ADD_LIBRARY( sll_fcisl STATIC sll_m_fcisl.F90 )
  TARGET_LINK_LIBRARIES( sll_fcisl sll_advection sll_utilities )

  ADD_LIBRARY( sll_derivative STATIC sll_m_derivative_2d_oblic.F90 )
  TARGET_LINK_LIBRARIES( sll_derivative 
    sll_advection
    sll_fcisl
    sll_utilities )

  ADD_LIBRARY( sll_fcisl_toroidal STATIC sll_m_fcisl_toroidal.F90 )
  TARGET_LINK_LIBRARIES( sll_fcisl_toroidal
    sll_advection
    sll_interpolators
    sll_utilities
    sll_file_io
    sll_parallel_array_utilities )


  # Tests
  IF(BUILD_TESTING)

    ADD_SUBDIRECTORY(testing)
  
    ADD_EXECUTABLE( test_aligned_derivative_2d aligned_derivative_2d.F90 )
    TARGET_LINK_LIBRARIES( test_aligned_derivative_2d 
      sll_derivative
      sll_fcisl
      sll_meshes 
      sll_advection
      sll_constants )


  ENDIF(BUILD_TESTING)

ENDIF(HDF5_ENABLED)
ENDIF(NOT STDF95)
