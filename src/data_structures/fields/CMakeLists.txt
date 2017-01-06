ADD_LIBRARY(sll_fields STATIC 
            sll_m_array_plotting.F90
            sll_m_scalar_field_1d_base.F90
            sll_m_scalar_field_1d.F90
            sll_m_scalar_field_2d_base.F90
            sll_m_scalar_field_2d.F90
            sll_m_dg_fields.F90)

TARGET_INCLUDE_DIRECTORIES(sll_fields BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})

TARGET_LINK_LIBRARIES(sll_fields
                      sll_coordinate_transformations
                      sll_meshes
                      sll_interpolators
                      sll_file_io
                      sll_utilities
                      sll_quadrature)

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


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
