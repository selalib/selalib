add_executable(test_hermite_interpolation_1d test_hermite_interpolation_1d.F90)
target_link_libraries(test_hermite_interpolation_1d sll_hermite_interpolation)

add_executable(test_hermite_interpolation_2d test_hermite_interpolation_2d.F90)
target_link_libraries(test_hermite_interpolation_2d sll_hermite_interpolation)

add_executable(test_hermite_aligned_interpolation_2d
               test_hermite_aligned_interpolation_2d.F90)
target_link_libraries(test_hermite_aligned_interpolation_2d
                      sll_hermite_interpolation)

# Ctest
add_test(NAME hermite_interpolation_1d COMMAND test_hermite_interpolation_1d)
set_tests_properties(hermite_interpolation_1d PROPERTIES PASS_REGULAR_EXPRESSION
                                                         "PASSED")

add_test(NAME hermite_interpolation_2d COMMAND test_hermite_interpolation_2d)
set_tests_properties(hermite_interpolation_2d PROPERTIES PASS_REGULAR_EXPRESSION
                                                         "PASSED")

add_test(NAME hermite_aligned_interpolation_2d
         COMMAND test_hermite_aligned_interpolation_2d)
set_tests_properties(hermite_aligned_interpolation_2d
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

# 2nd executable
add_executable(test_hex_hermite test_hex_hermite.F90)
target_link_libraries(test_hex_hermite sll_hermite_interpolation sll_meshes)

# 2nd bis executable ( circular advection )
add_executable(test_hex_hermite_circ test_hex_hermite_circ.F90)
target_link_libraries(test_hex_hermite_circ sll_hermite_interpolation
                      sll_meshes)

# 3rd executable
add_executable(test_hctc test_hctc.F90)
target_link_libraries(test_hctc sll_memory)

# 4th executable
add_executable(test_dimitrov test_dimitrov.F90)
target_link_libraries(test_dimitrov sll_memory)
