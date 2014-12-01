   ADD_LIBRARY(sll_distribution_function_multipatch
     STATIC sll_distribution_function_4d_multipatch.F90)

   TARGET_LINK_LIBRARIES(sll_distribution_function_multipatch
     sll_scalar_field_multipatch
     sll_coordinate_transformation_multipatch
     sll_remap
     sll_parallel_array_initializers
     sll_timer
     sll_logical_meshes_multipatch
     sll_memory)

   ADD_EXECUTABLE(test_distribution_function_4d_multipatch unit_test.F90)
   TARGET_LINK_LIBRARIES(test_distribution_function_4d_multipatch
     sll_distribution_function_multipatch
     sll_scalar_field_multipatch
     sll_common_array_initializers
     sll_timer
     sll_file_io
     sll_memory)
