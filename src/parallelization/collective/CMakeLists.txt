ADD_LIBRARY( sll_collective STATIC sll_m_collective.F90 )
TARGET_LINK_LIBRARIES( sll_collective
  sll_mpi
  sll_assert
  sll_memory
  sll_working_precision )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
