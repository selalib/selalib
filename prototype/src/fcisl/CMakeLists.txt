IF( NOT STDF95)
IF(HDF5_ENABLED)
#hdf5 is parallel and inside aligned_translation_2d.F90
#hdf5 is used
#so the compiler should be parallel
#should be changed
 
   IF(MPI_MODULE_ENABLED)
   
      SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

   ENDIF()

  ADD_LIBRARY (sll_fcisl STATIC
		fcisl.F90)

  TARGET_LINK_LIBRARIES(sll_fcisl sll_utilities sll_advection)

  ADD_EXECUTABLE(test_fcisl unit_test.F90)
  ADD_EXECUTABLE(test_aligned_derivative_2d aligned_derivative_2d.F90)

  TARGET_LINK_LIBRARIES(test_fcisl 
    sll_fcisl
    sll_logical_meshes 
    sll_advection
    sll_constants)

  TARGET_LINK_LIBRARIES(test_aligned_derivative_2d 
    sll_fcisl
    sll_logical_meshes 
    sll_advection
    sll_constants)

ENDIF(HDF5_ENABLED)

ENDIF(NOT STDF95)
