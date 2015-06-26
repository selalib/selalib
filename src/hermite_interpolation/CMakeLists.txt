ADD_LIBRARY ( sll_hermite_interpolation STATIC
  sll_hermite_interpolation_1d.F90 sll_hermite_interpolation_2d.F90 )
TARGET_LINK_LIBRARIES( sll_hermite_interpolation
  sll_utilities sll_boundary_condition_descriptors )

ADD_LIBRARY ( sll_hermite_aligned_interpolation STATIC
  sll_hermite_aligned_interpolation_2d.F90 )
TARGET_LINK_LIBRARIES( sll_hermite_aligned_interpolation
  sll_hermite_interpolation )


ADD_LIBRARY( hex_hermite STATIC sll_hex_interpolation_hermite.F90 )
TARGET_LINK_LIBRARIES( hex_hermite sll_hermite_interpolation sll_hex_meshes )


ADD_LIBRARY ( euler_2d_hex STATIC euler_2d_hex.F90)
TARGET_LINK_LIBRARIES( euler_2d_hex sll_utilities)

# Tests
IF(BUILD_TESTING)

  # 1st executable
  ADD_EXECUTABLE( test_hermite_interpolation unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_hermite_interpolation
    sll_hermite_interpolation
    sll_constants )

  # executable
  ADD_EXECUTABLE( test_hermite_aligned_interpolation test_hermite_aligned.F90 )
  TARGET_LINK_LIBRARIES( test_hermite_aligned_interpolation
    sll_hermite_aligned_interpolation
    sll_constants )


  # 2nd executable
  ADD_EXECUTABLE( test_hex_hermite test_hex_hermite.F90 )
  TARGET_LINK_LIBRARIES( test_hex_hermite
    sll_hermite_interpolation 
    hex_hermite 
    sll_hex_meshes 
    sll_constants )
  # 2nd bis executable ( circular advection )
  ADD_EXECUTABLE( test_hex_hermite_circ test_hex_hermite_circ.F90 )
  TARGET_LINK_LIBRARIES( test_hex_hermite_circ
    sll_hermite_interpolation 
    hex_hermite 
    sll_hex_meshes 
    sll_constants )

  # 3rd executable
  ADD_EXECUTABLE( test_hctc test_hctc.F90 )
  TARGET_LINK_LIBRARIES( test_hctc sll_constants sll_memory )

  # 4th executable
  ADD_EXECUTABLE( test_dimitrov test_dimitrov.F90 )
  TARGET_LINK_LIBRARIES( test_dimitrov sll_constants sll_memory )

  # 5th executable
  ADD_EXECUTABLE( test_hex_center_guide 
    test_hex_center_guide.F90 
    euler_2d_hex.F90
    ../utilities/pivotbande.F90 )
  TARGET_LINK_LIBRARIES( test_hex_center_guide 
    sll_hermite_interpolation 
    hex_hermite 
    sll_hex_meshes 
    hex_poisson
    sll_constants )



  
  # Ctest for executable 1
  ADD_TEST( NAME hermite_interpolation COMMAND test_hermite_interpolation )
  SET_TESTS_PROPERTIES( hermite_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  # Ctest for executable 1
  ADD_TEST( NAME hermite_aligned_interpolation COMMAND test_hermite_aligned_interpolation )
  SET_TESTS_PROPERTIES( hermite_aligned_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


ENDIF(BUILD_TESTING)
