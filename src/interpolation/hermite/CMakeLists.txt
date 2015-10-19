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


ENDIF(BUILD_TESTING)
