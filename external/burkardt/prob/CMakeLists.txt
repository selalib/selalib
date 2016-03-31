add_library( prob STATIC sll_m_prob.F90 )
install( TARGETS prob ARCHIVE DESTINATION lib )

# # Test
# if( BUILD_TESTING )
#   add_executable( prob_prb prob_prb.F90 )
#   target_link_libraries( prob_prb sll_prob )
# endif()
