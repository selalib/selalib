ADD_LIBRARY( sll_parallel_utilities STATIC buffer_loader_utilities.F90 )
TARGET_LINK_LIBRARIES( sll_parallel_utilities sll_remap sll_collective ) 

# Testing
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_buffer_loader_utilities
    unit_test_buffer_loader_utilities.F90 )
  TARGET_LINK_LIBRARIES( test_buffer_loader_utilities
    sll_parallel_utilities )
ENDIF(BUILD_TESTING)
