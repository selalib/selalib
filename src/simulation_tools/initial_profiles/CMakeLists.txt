ADD_LIBRARY(sll_initial_profiles STATIC 
            sll_m_common_array_initializers.F90 )
ADD_DEPENDENCIES(sll_initial_profiles sll_constants)

TARGET_LINK_LIBRARIES(sll_initial_profiles
                      sll_assert)
