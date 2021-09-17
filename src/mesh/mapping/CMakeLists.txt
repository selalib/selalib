IF(HDF5_IS_PARALLEL AND HDF5_ENABLED)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL AND HDF5_ENABLED)

ADD_LIBRARY(sll_mapping STATIC
  sll_m_mapping_2d.F90
  sll_m_mapping_3d.F90
  sll_m_3d_coordinate_transformations.F90)
  #sll_m_extraction_operator_0.F90
  #sll_m_extraction_operator_11.F90
  #sll_m_extraction_operator_12.F90
  #sll_m_extraction_operator_transposed_0.F90
  #sll_m_extraction_operator_transposed_11.F90
  #sll_m_extraction_operator_transposed_12.F90

TARGET_LINK_LIBRARIES( sll_mapping 
                       #sll_meshes
		       sll_splines
                       sll_linear_solvers)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
