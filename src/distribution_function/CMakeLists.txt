SET( DISTRIBUTION_SOURCES_LIBRARY 
               sll_scalar_field_1d_old.F90
               sll_scalar_field_2d_base_old.F90
               sll_scalar_field_2d_old.F90
               sll_fdistribu4d_DK.F90
	         initial_functions.F90
               landau_2d_initializer.F90
               tsi_2d_initializer.F90
               sll_distribution_function.F90
               scalar_field_initializer_base.F90)

IF(MPI_MODULE_ENABLED)

  ADD_LIBRARY(sll_distribution_function STATIC 
                 ${DISTRIBUTION_SOURCES_LIBRARY}
                 distribution_function_initializer_4d.F90
                 distribution_function_initializer_6d.F90)
  
  TARGET_LINK_LIBRARIES(sll_distribution_function sll_coordinate_transformations
                                                  sll_meshes
                                                  sll_file_io
                                                  sll_collective
                                                  sll_remap )
  # Ctest
  IF(BUILD_TESTING)
    ADD_EXECUTABLE(test_distribution_function unit_test.F90)
    TARGET_LINK_LIBRARIES(test_distribution_function sll_distribution_function)
    ADD_TEST( NAME distribution_function COMMAND test_distribution_function )
  ENDIF(BUILD_TESTING)

ELSE()

  ADD_LIBRARY(sll_distribution_function STATIC ${DISTRIBUTION_SOURCES_LIBRARY})
  
  TARGET_LINK_LIBRARIES(sll_distribution_function sll_coordinate_transformations
                                                  sll_meshes
                                                  sll_file_io)

ENDIF(MPI_MODULE_ENABLED)

