ADD_LIBRARY(sll_parallel_array_initializers STATIC 
            sll_parallel_array_initializer_module.F90
            sll_common_array_initializers_module.F90 )


TARGET_LINK_LIBRARIES(sll_parallel_array_initializers 
                      sll_remap
                      sll_coordinate_transformations
                      sll_logical_meshes)

ADD_LIBRARY(sll_common_array_initializers STATIC 
            sll_common_array_initializers_module.F90 )

TARGET_LINK_LIBRARIES(sll_common_array_initializers sll_assertion)

ADD_EXECUTABLE(test_parallel_array_initializer unit_test_4d.F90)

TARGET_LINK_LIBRARIES(test_parallel_array_initializer 
                         sll_parallel_array_initializers )
