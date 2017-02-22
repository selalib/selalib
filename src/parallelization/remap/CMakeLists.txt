ADD_LIBRARY( sll_remap STATIC sll_m_remapper.F90 )
TARGET_LINK_LIBRARIES( sll_remap sll_collective sll_utilities ) 


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
