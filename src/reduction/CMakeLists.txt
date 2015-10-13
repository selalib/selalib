ADD_LIBRARY( sll_reduction STATIC sll_reduction.F90 )
TARGET_LINK_LIBRARIES( sll_reduction sll_meshes ) 

# Ctest
IF(BUILD_TESTING)

  ADD_SUBDIRECTORY(testing)

ENDIF(BUILD_TESTING)
