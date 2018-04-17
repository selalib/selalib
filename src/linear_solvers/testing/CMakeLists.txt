ADD_EXECUTABLE( test_tridiagonal test_tridiagonal.F90 )
TARGET_LINK_LIBRARIES( test_tridiagonal sll_linear_solvers )
ADD_TEST( NAME tridiagonal COMMAND test_tridiagonal )

ADD_EXECUTABLE( test_toeplitz_penta_diagonal test_toeplitz_penta_diagonal.F90 )
TARGET_LINK_LIBRARIES( test_toeplitz_penta_diagonal sll_linear_solvers )
ADD_TEST( NAME toeplitz_penta_diagonal COMMAND test_toeplitz_penta_diagonal )
SET_TESTS_PROPERTIES( toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( example_tridiagonal example_tridiagonal.F90 )
TARGET_LINK_LIBRARIES( example_tridiagonal sll_linear_solvers )

#-------------------------------------------------------------------------------
# New tests
#-------------------------------------------------------------------------------

# conjugate gradient solver
add_executable( test_conjugate_gradient test_conjugate_gradient.F90 )
target_link_libraries( test_conjugate_gradient sll_vector_space sll_linear_solvers )
add_test( NAME conjugate_gradient COMMAND test_conjugate_gradient )
set_tests_properties( conjugate_gradient PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

# dot product (matrix-vector) in stencil format
add_executable( test_dot_product_stencil test_dot_product_stencil.F90 )
target_link_libraries( test_dot_product_stencil sll_vector_space sll_linear_solvers )
add_test( NAME dot_product_stencil COMMAND test_dot_product_stencil )
set_tests_properties( dot_product_stencil PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
