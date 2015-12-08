ADD_LIBRARY(sll_common_array_initializers STATIC 
            sll_m_common_array_initializers.F90 )
ADD_DEPENDENCIES(sll_common_array_initializers sll_constants)

TARGET_LINK_LIBRARIES(sll_common_array_initializers
                      sll_assert)
