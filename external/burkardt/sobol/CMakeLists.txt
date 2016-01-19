add_library( sobol STATIC sll_m_sobol.F90 )
install( TARGETS sobol ARCHIVE DESTINATION lib )

# # Tests
# if( BUILD_TESTING )
#   add_executable( test_sobol_prb sobol_prb.f90 )
#   target_link_libraries( test_sobol_prb sobol )
# endif()
