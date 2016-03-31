ADD_EXECUTABLE( test_hermite_interpolation_1d test_hermite_interpolation_1d.F90 )
TARGET_LINK_LIBRARIES( test_hermite_interpolation_1d
  sll_hermite_interpolation)

ADD_EXECUTABLE( test_hermite_interpolation_2d test_hermite_interpolation_2d.F90 )
TARGET_LINK_LIBRARIES( test_hermite_interpolation_2d
  sll_hermite_interpolation)

ADD_EXECUTABLE( test_hermite_aligned_interpolation_2d 
  test_hermite_aligned_interpolation_2d.F90 )
TARGET_LINK_LIBRARIES( test_hermite_aligned_interpolation_2d
  sll_hermite_interpolation)

# Ctest
ADD_TEST( NAME hermite_interpolation_1d COMMAND test_hermite_interpolation_1d )
SET_TESTS_PROPERTIES( hermite_interpolation_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_TEST( NAME hermite_interpolation_2d COMMAND test_hermite_interpolation_2d )
SET_TESTS_PROPERTIES( hermite_interpolation_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_TEST( NAME hermite_aligned_interpolation_2d COMMAND test_hermite_aligned_interpolation_2d )
SET_TESTS_PROPERTIES( hermite_aligned_interpolation_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

# 2nd executable
ADD_EXECUTABLE( test_hex_hermite test_hex_hermite.F90 )
TARGET_LINK_LIBRARIES( test_hex_hermite
  sll_hermite_interpolation
  sll_meshes)

# 2nd bis executable ( circular advection )
ADD_EXECUTABLE( test_hex_hermite_circ test_hex_hermite_circ.F90 )
TARGET_LINK_LIBRARIES( test_hex_hermite_circ
  sll_hermite_interpolation
  sll_meshes)

# 3rd executable
ADD_EXECUTABLE( test_hctc test_hctc.F90 )
TARGET_LINK_LIBRARIES( test_hctc  sll_memory )

# 4th executable
ADD_EXECUTABLE( test_dimitrov test_dimitrov.F90 )
TARGET_LINK_LIBRARIES( test_dimitrov  sll_memory )

# 5th executable
ADD_EXECUTABLE( test_hex_center_guide
  test_hex_center_guide.F90)

TARGET_LINK_LIBRARIES( test_hex_center_guide
  sll_hermite_interpolation
  sll_meshes
  sll_poisson_solvers
  sll_utilities )
