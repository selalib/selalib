# Test vector space wrapper for 1D arrays
add_executable( test_vector_space_real_array_1d test_vector_space_real_array_1d.F90 )
target_link_libraries( test_vector_space_real_array_1d sll_vector_space )
add_test( NAME vector_space_real_array_1d COMMAND test_vector_space_real_array_1d )
set_tests_properties( vector_space_real_array_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

# Test vector space wrapper for 2D arrays
add_executable( test_vector_space_real_array_2d test_vector_space_real_array_2d.F90 )
target_link_libraries( test_vector_space_real_array_2d sll_vector_space )
add_test( NAME vector_space_real_array_2d COMMAND test_vector_space_real_array_2d )
set_tests_properties( vector_space_real_array_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


# Test vector space wrapper for 3D arrays
add_executable( test_vector_space_real_array_3d test_vector_space_real_array_3d.F90 )
target_link_libraries( test_vector_space_real_array_3d sll_vector_space )
add_test( NAME vector_space_real_array_3d COMMAND test_vector_space_real_array_3d )
set_tests_properties( vector_space_real_array_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
