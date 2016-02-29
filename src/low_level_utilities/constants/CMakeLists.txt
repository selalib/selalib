ADD_LIBRARY( sll_constants STATIC sll_m_constants.F90 )
TARGET_LINK_LIBRARIES( sll_constants sll_working_precision)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
