ADD_LIBRARY( sll_constants OBJECT sll_m_constants.F90 )
ADD_DEPENDENCIES( sll_constants sll_working_precision)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
