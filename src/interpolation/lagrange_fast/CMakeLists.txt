ADD_LIBRARY( sll_lagrange_fast STATIC
             sll_m_lagrange_fast.F90 )

#TARGET_LINK_LIBRARIES( sll_lagrange_interpolation sll_utilities sll_boundary_condition_descriptors)

# Ctest
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
