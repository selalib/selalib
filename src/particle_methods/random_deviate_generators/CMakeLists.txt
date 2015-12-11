add_library( sll_random_deviate_generators STATIC sll_m_gaussian.F90
						  sll_m_hammersley.F90 
						  #sll_m_marsaglia.F90 
)
target_link_libraries( sll_random_deviate_generators
sll_assert sll_memory )

ADD_SUBDIRECTORY(moment_matching)


# Test
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
