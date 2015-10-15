ADD_LIBRARY( sll_collective STATIC sll_collective.F90 )
TARGET_LINK_LIBRARIES( sll_collective sll_utilities )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
