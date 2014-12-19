ADD_LIBRARY(sll_fields STATIC 
            sll_scalar_field_1d_base.F90
            sll_scalar_field_1d.F90
            sll_scalar_field_2d_base.F90
            sll_scalar_field_2d.F90)		

TARGET_LINK_LIBRARIES(sll_fields
                      sll_coordinate_transformations
                      sll_meshes
                      sll_interpolators
                      sll_file_io
                      sll_utilities)

#IF (MPI_ENABLED)
#   SET(SLL_SCALAR_FIELD_SRCS ${SLL_SCALAR_FIELD_SRCS}
#                             distribution_function_initializer_4d.F90
#                             distribution_function_initializer_6d.F90)
#ENDIF (MPI_ENABLED)

#ADD_LIBRARY(sll_scalar_field STATIC ${SLL_SCALAR_FIELD_SRCS} )

#TARGET_LINK_LIBRARIES( sll_scalar_field sll_coordinate_transformations
#                                        sll_utilities)
#
#ADD_EXECUTABLE(test_scalar_field unit_test.F90)
#
#TARGET_LINK_LIBRARIES(test_scalar_field sll_scalar_field)

ADD_LIBRARY(sll_dg_fields STATIC sll_dg_fields_2d.F90 
                                 init_function_module.F90)		

TARGET_LINK_LIBRARIES(sll_dg_fields sll_coordinate_transformations 
                                    sll_integration)

ADD_LIBRARY(sll_fem_fields STATIC sll_fem_fields_2d.F90 
                                 init_function_module.F90)		

TARGET_LINK_LIBRARIES(sll_fem_fields sll_coordinate_transformations 
                                    sll_integration)


# Tests
IF(BUILD_TESTING)

  # Executables not added to ctest
  ADD_EXECUTABLE(test_dg_fields unit_test_dg_fields.F90)
  TARGET_LINK_LIBRARIES(test_dg_fields sll_dg_fields)
  ADD_EXECUTABLE(test_fem_fields unit_test_fem_fields.F90)
  TARGET_LINK_LIBRARIES(test_fem_fields sll_fem_fields)

  # Ctest 1D
  ADD_EXECUTABLE( test_scalar_field_1d unit_test_1d.F90  )
  TARGET_LINK_LIBRARIES( test_scalar_field_1d sll_fields )
  ADD_TEST( NAME fields_1d COMMAND test_scalar_field_1d )
  SET_TESTS_PROPERTIES( fields_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  # Ctest 2D
  ADD_EXECUTABLE( test_scalar_field_2d unit_test_2d.F90 helper_functions.F90 )
  TARGET_LINK_LIBRARIES( test_scalar_field_2d sll_fields )
  ADD_TEST( NAME fields_2d COMMAND test_scalar_field_2d )
  SET_TESTS_PROPERTIES( fields_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(BUILD_TESTING)
