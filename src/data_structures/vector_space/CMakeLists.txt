add_library( sll_vector_space STATIC
             sll_m_vector_space_base.F90
             sll_m_vector_space_real_arrays.F90 )

target_link_libraries( sll_vector_space
                       sll_errors )
# Tests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )
