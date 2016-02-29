ADD_LIBRARY( sll_collective STATIC sll_m_collective.F90 )
TARGET_LINK_LIBRARIES( sll_collective
  sll_assert
  sll_memory
  sll_mpi )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
