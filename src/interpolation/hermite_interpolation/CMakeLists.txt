ADD_LIBRARY ( sll_hermite_interpolation STATIC
              sll_m_hermite_interpolation_1d.F90
              sll_m_hermite_interpolation_2d.F90
              sll_m_interpolation_hex_hermite.F90
              sll_m_hermite_aligned_interpolation_2d.F90 )

TARGET_LINK_LIBRARIES( sll_hermite_interpolation
                         sll_meshes
                         sll_utilities)

# Ctests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
