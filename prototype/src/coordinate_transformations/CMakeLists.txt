ADD_LIBRARY(sll_coordinate_transformations STATIC
            sll_coordinate_transformation_2d_base.F90
            sll_coordinate_transformations_2d.F90
	    sll_coordinate_transformations_2d_nurbs.F90
            sll_common_coordinate_transformations.F90)

TARGET_LINK_LIBRARIES(sll_coordinate_transformations sll_interpolators
	                                               sll_file_io
	                                               sll_logical_meshes)

IF(HDF5_IS_PARALLEL)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

ADD_EXECUTABLE(test_coordinate_transformations_2d unit_test_2d.F90)

TARGET_LINK_LIBRARIES(test_coordinate_transformations_2d 
                         sll_coordinate_transformations)

