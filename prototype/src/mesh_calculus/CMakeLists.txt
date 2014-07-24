IF(HDF5_IS_PARALLEL AND HDF5_ENABLED)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL AND HDF5_ENABLED)

ADD_LIBRARY(sll_mesh_calculus STATIC sll_mesh_calculus_module.F90 )

TARGET_LINK_LIBRARIES(sll_mesh_calculus sll_coordinate_transformations
                                        sll_logical_meshes
                                        sll_integration
                                        sll_assertion
                                        sll_working_precision)

ADD_EXECUTABLE( test_mesh_calculus_2d unit_test_2d.F90 )

TARGET_LINK_LIBRARIES( test_mesh_calculus_2d sll_mesh_calculus sll_file_io sll_memory ${HDF5_LIBRARIES} )


