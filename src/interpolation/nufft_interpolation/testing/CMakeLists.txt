add_executable(test_nufft_interpolation test_nufft_interpolation.F90)
target_link_libraries(test_nufft_interpolation sll_nufft_interpolation
                      sll_file_io)
add_test(NAME nufft_interpolation COMMAND test_nufft_interpolation)
set_tests_properties(nufft_interpolation PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")
