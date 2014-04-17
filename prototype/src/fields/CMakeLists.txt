ADD_LIBRARY(sll_scalar_field_2d_alternative STATIC 
            sll_scalar_field_2d_base.F90
            sll_scalar_field_2d_alternative.F90)		

TARGET_LINK_LIBRARIES(sll_scalar_field_2d_alternative
                      sll_coordinate_transformations
                      sll_logical_meshes
                      sll_interpolators
                      sll_file_io
                      sll_utilities)

ADD_EXECUTABLE(test_scalar_field_alternative 
  unit_test_alternative.F90 
  helper_functions_unit_test_alternative.F90)

TARGET_LINK_LIBRARIES(test_scalar_field_alternative sll_scalar_field_2d_alternative)

SET(SLL_SCALAR_FIELD_SRCS sll_scalar_field_1d.F90
                          scalar_field_initializer_base.F90
                          sll_scalar_field_2d_base.F90
                          sll_scalar_field_2d.F90 
                          sll_scalar_field_2d_alternative.F90
                          scalar_field_initializer_base.F90
                          landau_2d_initializer.F90
                          gaussian_2d_initializer.F90
                          tsi_2d_initializer.F90)


ADD_LIBRARY(sll_scalar_field_1d_alternative STATIC 
               sll_scalar_field_1d_base_alternative.F90
               sll_scalar_field_1d_alternative.F90)		

TARGET_LINK_LIBRARIES(sll_scalar_field_1d_alternative
                      sll_logical_meshes
                      sll_interpolators
                      sll_file_io
                      sll_utilities)

ADD_EXECUTABLE(test_scalar_fields_1d_alternative unit_test_1d_alternative.F90)

TARGET_LINK_LIBRARIES(test_scalar_fields_1d_alternative sll_scalar_field_1d_alternative)

IF (MPI_MODULE_ENABLED)
   SET(SLL_SCALAR_FIELD_SRCS ${SLL_SCALAR_FIELD_SRCS}
                             distribution_function_initializer_4d.F90
                             distribution_function_initializer_6d.F90)
ENDIF (MPI_MODULE_ENABLED)

ADD_LIBRARY(sll_scalar_field STATIC ${SLL_SCALAR_FIELD_SRCS} )

TARGET_LINK_LIBRARIES( sll_scalar_field sll_mapped_meshes
                                        sll_coordinate_transformations
                                        sll_utilities)

ADD_EXECUTABLE(test_scalar_field unit_test.F90)

TARGET_LINK_LIBRARIES(test_scalar_field sll_scalar_field)

ADD_LIBRARY(sll_dg_fields STATIC sll_dg_fields_2d.F90)		

TARGET_LINK_LIBRARIES(sll_dg_fields sll_integration)
