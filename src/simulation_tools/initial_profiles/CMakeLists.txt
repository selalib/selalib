add_library(sll_initial_profiles STATIC sll_m_common_array_initializers.F90)

target_link_libraries(sll_initial_profiles sll_constants sll_assert
                      sll_random_deviate_generators sll_errors)
