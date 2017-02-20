add_executable( test_random_deviate_generators test_random_deviate_generators.F90 )
target_link_libraries( test_random_deviate_generators
  sll_random_deviate_generators
  sll_assert
  sll_memory)


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

FIND_PACKAGE(GSL)

IF(FGSL_FOUND AND NOT PGI)
  ADD_EXECUTABLE( test_fgsl_rng mod_unit.f90 rng.f90)
  TARGET_INCLUDE_DIRECTORIES( test_fgsl_rng PRIVATE ${FGSL_INCLUDES})
  TARGET_LINK_LIBRARIES( test_fgsl_rng ${FGSL_LIBRARIES} )
ENDIF()
