SET ( DG_INTERPOLATION_SRC_FILES
  sll_m_dg_interpolator_1d.F90
  )

SET (DG_INTERPOLATION_DEPS
  sll_quadrature
  sll_constants
  sll_interpolators )

ADD_LIBRARY( sll_dg_interpolation STATIC ${DG_INTERPOLATION_SRC_FILES})

TARGET_LINK_LIBRARIES( sll_dg_interpolation ${DG_INTERPOLATION_DEPS} )


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
