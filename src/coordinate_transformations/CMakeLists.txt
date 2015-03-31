IF(HDF5_ENABLED AND HDF5_IS_PARALLEL)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL)

ADD_LIBRARY(sll_coordinate_transformations STATIC
            sll_coordinate_transformation_2d_base.F90
            sll_coordinate_transformations_2d.F90
            sll_coordinate_transformations_2d_nurbs.F90
            sll_common_coordinate_transformations.F90)

TARGET_LINK_LIBRARIES(sll_coordinate_transformations 
                      sll_interpolators sll_file_io sll_meshes)

IF(BUILD_TESTING)

  ADD_EXECUTABLE( test_coordinate_transformations_2d unit_test_2d.F90)

  TARGET_LINK_LIBRARIES( test_coordinate_transformations_2d 
    sll_coordinate_transformations )

  ADD_TEST( NAME coordinate_transformations COMMAND test_coordinate_transformations_2d )

  SET_TESTS_PROPERTIES( coordinate_transformations 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  SET( INPUT_FILES_DIR
    ${CMAKE_SOURCE_DIR}/CAID_coordinate_transformations/mesh_collela/ )

  FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )

  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ENDIF(BUILD_TESTING)

