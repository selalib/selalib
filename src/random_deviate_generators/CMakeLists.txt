add_library( sll_random_deviate_generators STATIC gaussian.F90 hammersley.F90 )
target_link_libraries( sll_random_deviate_generators
  sll_assert sll_memory sll_working_precision )

ADD_SUBDIRECTORY(sobol) #sobol sequence
ADD_SUBDIRECTORY(prob) #Probability density functions

add_library( sll_moment_matching STATIC moment_matching.F90)
target_link_libraries(sll_moment_matching sll_assert
					  sll_memory
					  sll_constants
					  sll_working_precision)	
# Test
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)