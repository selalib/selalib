IF(STDF95)
 ADD_LIBRARY(sll_mapped_meshes STATIC
	sll_mapped_meshes_1d.F90
	sll_mapped_meshes_2d.F90
	sll_mapped_meshes_2d_cartesian.F90
  )
ELSE()
 ADD_LIBRARY(sll_mapped_meshes STATIC
	sll_mapped_meshes_1d_base.F90
	sll_mapped_meshes_2d_base.F90
	sll_mapped_meshes_1d.F90
	sll_mapped_meshes_2d.F90
	sll_mapped_meshes_2d_cartesian.F90
  )
ENDIF()

TARGET_LINK_LIBRARIES(sll_mapped_meshes sll_interpolators 
                                        sll_coordinate_transformations
                                        sll_file_io)

ADD_EXECUTABLE(test_mapped_meshes_1d unit_test_1D.F90)

TARGET_LINK_LIBRARIES(test_mapped_meshes_1d sll_mapped_meshes)

IF(HDF5_IS_PARALLEL)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

ADD_EXECUTABLE(test_mapped_meshes_2d unit_test.F90)

TARGET_LINK_LIBRARIES(test_mapped_meshes_2d sll_mapped_meshes)

IF(NOT STDF95)

 ADD_EXECUTABLE(test_mapped_meshes_cartesian unit_test_cartesian.F90)
 TARGET_LINK_LIBRARIES(test_mapped_meshes_cartesian sll_mapped_meshes)

ENDIF()
