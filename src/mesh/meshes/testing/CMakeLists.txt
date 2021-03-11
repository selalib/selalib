# Cartesian mesh ..........
add_executable(test_cartesian_meshes test_cartesian_meshes.F90)
target_link_libraries(test_cartesian_meshes sll_meshes)
# Ctest:
add_test(NAME cartesian_meshes COMMAND test_cartesian_meshes)
set_tests_properties(cartesian_meshes PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED")
# ..........................

# Hexagonal mesh ..........
add_executable(test_hexagonal_meshes test_hexagonal_meshes.F90)
target_link_libraries(test_hexagonal_meshes sll_meshes)
# Ctest:
add_test(NAME hexagonal_meshes COMMAND test_hexagonal_meshes)
set_tests_properties(hexagonal_meshes PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED")
# ..........................

# Triangular mesh .........
add_executable(test_triangular_meshes test_triangular_meshes.F90)
target_link_libraries(test_triangular_meshes sll_meshes)
# Ctest:
set(REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/diode.maa)
add_test(NAME triangular_meshes COMMAND test_triangular_meshes
                                        ${REFERENCE_FILE})
set_tests_properties(triangular_meshes PROPERTIES PASS_REGULAR_EXPRESSION
                                                  "PASSED")
# ..........................
