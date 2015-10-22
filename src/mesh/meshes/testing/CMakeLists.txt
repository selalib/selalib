# Cartesian mesh ..........
ADD_EXECUTABLE( test_cartesian_meshes test_cartesian_meshes.F90 )
TARGET_LINK_LIBRARIES( test_cartesian_meshes sll_meshes )
# Ctest:
ADD_TEST( NAME cartesian_meshes COMMAND test_cartesian_meshes )
SET_TESTS_PROPERTIES( cartesian_meshes PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#..........................

# Hexagonal mesh ..........
ADD_EXECUTABLE( test_hexagonal_meshes test_hexagonal_meshes.F90 )
TARGET_LINK_LIBRARIES( test_hexagonal_meshes sll_meshes )
# Ctest:
ADD_TEST( NAME hexagonal_meshes COMMAND test_hexagonal_meshes )
SET_TESTS_PROPERTIES( hexagonal_meshes PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#..........................

# Triangular mesh .........
ADD_EXECUTABLE( test_triangular_meshes test_triangular_meshes.F90 )
TARGET_LINK_LIBRARIES( test_triangular_meshes sll_meshes )
# Ctest:
SET( REFERENCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/diode.maa )
ADD_TEST( NAME triangular_meshes COMMAND test_triangular_meshes ${REFERENCE_FILE})
SET_TESTS_PROPERTIES( triangular_meshes PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
#..........................
