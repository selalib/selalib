add_library( sll_vector_space STATIC
             sll_vector_space_base.F90
             sll_vector_space_real_arrays.F90 )

target_link_libraries( sll_vector_space
                       sll_working_precision
                       sll_errors )
# Tests
if( BUILD_TESTING )
  add_executable( test_vector_space test_vector_space.F90 )
  target_link_libraries( test_vector_space sll_vector_space )
endif( BUILD_TESTING )
