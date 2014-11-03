ADD_LIBRARY(sll_distribution_function STATIC 
               sll_scalar_field_2d_base_old.F90
               sll_scalar_field_2d_old.F90
               sll_fdistribu4d_DK.F90
	         initial_functions.F90
               landau_2d_initializer.F90
               tsi_2d_initializer.F90
               sll_distribution_function.F90
               distribution_function_initializer_4d.F90
               distribution_function_initializer_6d.F90
               scalar_field_initializer_base.F90)

TARGET_LINK_LIBRARIES(sll_distribution_function sll_coordinate_transformations
                                                 sll_meshes
                                                 sll_file_io)

ADD_EXECUTABLE(test_distribution_function unit_test.F90)

TARGET_LINK_LIBRARIES(test_distribution_function sll_distribution_function)
