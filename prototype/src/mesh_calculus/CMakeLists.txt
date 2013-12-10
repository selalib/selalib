IF(HDF5_IS_PARALLEL)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

ADD_LIBRARY(sll_mesh_calculus STATIC sll_mesh_calculus_module.F90 )

ADD_DEPENDENCIES(sll_mesh_calculus
  sll_coordinate_transformations
  sll_utilities
  sll_file_io
  sll_logical_meshes
  sll_integration
  sll_assertion
  sll_working_precision
  ${HDF5_LIBRARIES} )

ADD_EXECUTABLE( test_mesh_calculus_2d unit_test_2d.F90 )

TARGET_LINK_LIBRARIES( test_mesh_calculus_2d
  sll_mesh_calculus
  sll_coordinate_transformations
  sll_utilities
  sll_logical_meshes
  sll_file_io
  sll_integration
  sll_memory
  sll_assertion
  sll_working_precision
  ${HDF5_LIBRARIES} )


