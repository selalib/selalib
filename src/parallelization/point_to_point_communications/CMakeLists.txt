ADD_LIBRARY( sll_point_to_point_communications STATIC
  sll_m_point_to_point_comms.F90 )

TARGET_LINK_LIBRARIES( sll_point_to_point_communications
  sll_mpi sll_collective )


# Ctests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
