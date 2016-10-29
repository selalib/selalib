SET( CHARACTERISTICS_SRC_FILES 
  sll_m_characteristics_1d_explicit_euler.F90
  sll_m_characteristics_1d_explicit_euler_conservative.F90
  sll_m_characteristics_1d_trapezoid.F90
  sll_m_characteristics_1d_trapezoid_conservative.F90
  sll_m_characteristics_1d_base.F90
  sll_m_characteristics_2d_verlet.F90
  sll_m_characteristics_2d_explicit_euler.F90
  sll_m_characteristics_2d_explicit_euler_conservative.F90
  sll_m_characteristics_2d_base.F90
  sll_m_euler_2d_hex.F90 )

ADD_LIBRARY( sll_characteristics STATIC ${CHARACTERISTICS_SRC_FILES} )
TARGET_LINK_LIBRARIES( sll_characteristics sll_utilities sll_interpolators )

# Ctests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
