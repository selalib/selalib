ADD_LIBRARY ( sll_hermite_interpolation STATIC
  sll_m_hermite_interpolation_1d.F90 sll_m_hermite_interpolation_2d.F90 )
TARGET_LINK_LIBRARIES( sll_hermite_interpolation
  sll_utilities sll_boundary_condition_descriptors )

ADD_LIBRARY ( sll_hermite_aligned_interpolation STATIC
  sll_m_hermite_aligned_interpolation_2d.F90 )
TARGET_LINK_LIBRARIES( sll_hermite_aligned_interpolation
  sll_hermite_interpolation )


ADD_LIBRARY( hex_hermite STATIC sll_m_interpolation_hex_hermite.F90 )
TARGET_LINK_LIBRARIES( hex_hermite sll_hermite_interpolation sll_meshes )


ADD_LIBRARY ( euler_2d_hex STATIC sll_m_euler_2d_hex.F90)
TARGET_LINK_LIBRARIES( euler_2d_hex sll_utilities)

# Tests
IF(BUILD_TESTING)

  # Everything must go inside here!
  ADD_SUBDIRECTORY(testing)

  # 2nd executable
  ADD_EXECUTABLE( test_hex_hermite test_hex_hermite.F90 )
  TARGET_LINK_LIBRARIES( test_hex_hermite
    sll_hermite_interpolation
    hex_hermite
    sll_meshes
    sll_constants )

  # 2nd bis executable ( circular advection )
  ADD_EXECUTABLE( test_hex_hermite_circ test_hex_hermite_circ.F90 )
  TARGET_LINK_LIBRARIES( test_hex_hermite_circ
    sll_hermite_interpolation
    hex_hermite
    sll_meshes
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
    sll_m_euler_2d_hex.F90 )
  TARGET_LINK_LIBRARIES( test_hex_center_guide
    sll_hermite_interpolation
    hex_hermite
    sll_meshes
    hex_poisson
    sll_constants
    sll_utilities )

ENDIF(BUILD_TESTING)
