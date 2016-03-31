ADD_LIBRARY(sll_initial_profiles STATIC 
            sll_m_common_array_initializers.F90 )

TARGET_LINK_LIBRARIES(sll_initial_profiles
                      sll_constants
                      sll_assert)
