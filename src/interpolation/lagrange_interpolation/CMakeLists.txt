ADD_LIBRARY( sll_lagrange_interpolation STATIC
             sll_m_lagrange_interpolation_1d.F90)

TARGET_LINK_LIBRARIES( sll_lagrange_interpolation sll_utilities)

# Ctest
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
