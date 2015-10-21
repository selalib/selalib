ADD_LIBRARY( sll_parallel_utilities STATIC sll_m_buffer_loader_utilities.F90 )
TARGET_LINK_LIBRARIES( sll_parallel_utilities sll_remap sll_collective ) 
   

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)

