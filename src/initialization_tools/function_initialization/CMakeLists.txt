ADD_LIBRARY(sll_common_array_initializers STATIC 
            sll_m_common_array_initializers.F90 )

TARGET_LINK_LIBRARIES(sll_common_array_initializers
                      sll_assert)
