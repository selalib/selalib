ADD_LIBRARY( sll_decomposition STATIC sll_decomposition.F90 )
TARGET_LINK_LIBRARIES( sll_decomposition sll_collective )

# Ctests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)