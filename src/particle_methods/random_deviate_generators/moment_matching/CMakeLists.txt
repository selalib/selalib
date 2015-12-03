add_library( sll_moment_matching STATIC sll_m_moment_matching.F90 )
target_link_libraries(sll_moment_matching sll_assert
					  sll_memory)
# Test
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
 
