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

IF(BUILD_SPRNG)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  ADD_DEFINITIONS(-DPOINTER_SIZE=8 -DLONG64=long -DINTEGER_STAR_8 -DUSE_MPI)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/sprng/include)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/sprng/lib)
  ADD_EXECUTABLE(test_sprng test_sprng.F90)
  ADD_EXECUTABLE(test_sprng_mpi 2streamsf_mpi.F90)
  TARGET_LINK_LIBRARIES(test_sprng sprng stdc++ ${MPI_CXX_LIBRARIES})
  TARGET_LINK_LIBRARIES(test_sprng_mpi sprng stdc++ ${MPI_CXX_LIBRARIES})
ENDIF(BUILD_SPRNG)
