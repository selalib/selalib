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
  add_executable( test_random_deviate_generators unit_test.F90 )
  target_link_libraries( test_random_deviate_generators
    sll_random_deviate_generators
    sll_assert
    sll_memory
    sll_working_precision )
ENDIF(BUILD_TESTING)


INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/sprng/include)
ADD_EXECUTABLE(test_sprng 2streamsf_mpi.F90)
ADD_DEPENDENCIES(test_sprng sprng)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/sprng/lib)
TARGET_LINK_LIBRARIES(test_sprng sprng)
