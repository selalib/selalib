IF(HDF5_IS_PARALLEL AND HDF5_ENABLED)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL AND HDF5_ENABLED)

ADD_LIBRARY(sll_mesh_calculus STATIC sll_mesh_calculus_module.F90 )

TARGET_LINK_LIBRARIES(sll_mesh_calculus sll_coordinate_transformations
                                        sll_meshes
                                        sll_integration
                                        sll_assert
                                        sll_working_precision)

# Tests
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_mesh_calculus_2d unit_test_2d.F90 )
  TARGET_LINK_LIBRARIES( test_mesh_calculus_2d
    sll_mesh_calculus sll_file_io sll_memory ${HDF5_LIBRARIES} )
ENDIF(BUILD_TESTING)


SET(f90_sources poisson.F90 
                maillage.F90 solveurs_module.F90 
                errout.F90
                zone.F90 utlib.F90)

ADD_LIBRARY(picsou ${f90_sources})
ADD_EXECUTABLE(test_picsou picsou.F90)
TARGET_LINK_LIBRARIES(test_picsou picsou sll_meshes sll_file_io)
