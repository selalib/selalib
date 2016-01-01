ADD_LIBRARY( sll_timer STATIC sll_m_timer.F90 )
ADD_DEPENDENCIES(sll_timer sll_working_precision)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
