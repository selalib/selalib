add_executable( test_random_deviate_generators test_random_deviate_generators.F90 )
target_link_libraries( test_random_deviate_generators
  sll_random_deviate_generators
  sll_assert
  sll_memory)



FIND_PACKAGE(GSL)

IF(FGSL_FOUND AND NOT PGI)
  ADD_EXECUTABLE( test_fgsl_rng mod_unit.f90 rng.f90)
  TARGET_INCLUDE_DIRECTORIES( test_fgsl_rng PRIVATE ${FGSL_INCLUDES})
  TARGET_LINK_LIBRARIES( test_fgsl_rng ${FGSL_LIBRARIES} )
ENDIF()
