IF(HDF5_IS_PARALLEL)
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)
ADD_SUBDIRECTORY(CAID_coordinate_transformations)
ADD_SUBDIRECTORY(CAID_coordinate_transformations/square_4p_n10)

ADD_EXECUTABLE( test_distribution_function_4d_multipatch
                test_distribution_function_4d_multipatch.F90 )

TARGET_LINK_LIBRARIES( test_distribution_function_4d_multipatch sll_multipatch
  sll_common_array_initializers
  sll_timer
  sll_file_io)

IF(PYTHON3_FOUND)
  SET( ARGS " " )
  ADD_MPI_TEST( distribution_function_4d_multipatch 
     test_distribution_function_4d_multipatch ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( distribution_function_4d_multipatch 
     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(PYTHON3_FOUND)

SET( INPUT_FILES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/square_4p_n10/ )
FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ADD_EXECUTABLE( test_meshes_multipatch test_meshes_multipatch.F90)
TARGET_LINK_LIBRARIES( test_meshes_multipatch sll_multipatch )
ADD_TEST( NAME meshes_multipatch COMMAND test_meshes_multipatch )
SET_TESTS_PROPERTIES( meshes_multipatch PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_coordinate_transformation_multipatch_2d 
                test_coordinate_transformation_multipatch_2d.F90 )

TARGET_LINK_LIBRARIES(test_coordinate_transformation_multipatch_2d
   sll_multipatch sll_file_io)

IF(PYTHON3_FOUND)
  ADD_TEST( NAME coordinate_transformation_multipatch_2d
  COMMAND test_coordinate_transformation_multipatch_2d )
  SET_TESTS_PROPERTIES( coordinate_transformation_multipatch_2d 
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(PYTHON3_FOUND)

SET( INPUT_FILES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/circle_mp5_pts12/ )
FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ADD_EXECUTABLE( test_scalar_field_multipatch_2d
                test_scalar_field_multipatch_2d.F90 )
TARGET_LINK_LIBRARIES( test_scalar_field_multipatch_2d sll_multipatch )

IF(PYTHON3_FOUND)
  ADD_TEST( NAME scalar_field_multipatch_2d COMMAND test_scalar_field_multipatch_2d )
  SET_TESTS_PROPERTIES( scalar_field_multipatch_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(PYTHON3_FOUND)

SET( INPUT_FILES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/square_4p_n10/ )
FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

ADD_EXECUTABLE( test_general_coordinate_elliptic_solver_multipatch 
                m_multipatch_helper_functions.F90 
                test_general_coordinate_elliptic_solver_multipatch.F90)

TARGET_LINK_LIBRARIES( test_general_coordinate_elliptic_solver_multipatch
  sll_multipatch)

IF(SUITESPARSE_ENABLED)
ELSE()
IF(PYTHON3_FOUND)
  ADD_TEST( NAME general_coordinate_elliptic_solver_multipatch 
  COMMAND test_general_coordinate_elliptic_solver_multipatch )
  SET_TESTS_PROPERTIES( general_coordinate_elliptic_solver_multipatch
  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(PYTHON3_FOUND)
ENDIF(SUITESPARSE_ENABLED)

SET( INPUT_FILES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/circle_mp5_pts12/ )
FILE( GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml" )
FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
