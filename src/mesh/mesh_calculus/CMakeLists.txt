IF(HDF5_IS_PARALLEL AND HDF5_ENABLED)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL AND HDF5_ENABLED)

ADD_LIBRARY(sll_mesh_calculus STATIC sll_m_mesh_calculus_2d.F90 )

TARGET_LINK_LIBRARIES(sll_mesh_calculus sll_coordinate_transformations
                                        sll_meshes
                                        sll_quadrature
                                        sll_assert)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
