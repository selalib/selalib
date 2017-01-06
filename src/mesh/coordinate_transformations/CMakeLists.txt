IF(HDF5_ENABLED AND HDF5_IS_PARALLEL)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL)

ADD_LIBRARY(sll_coordinate_transformations STATIC
            sll_m_coordinate_transformation_2d_base.F90
            sll_m_coordinate_transformations_2d.F90
            sll_m_common_coordinate_transformations.F90)

TARGET_INCLUDE_DIRECTORIES(sll_coordinate_transformations BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})
TARGET_LINK_LIBRARIES(sll_coordinate_transformations 
                      sll_interpolators sll_file_io sll_meshes)
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
