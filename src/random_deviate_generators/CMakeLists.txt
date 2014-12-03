add_library( sll_random_deviate_generators STATIC gaussian.F90 hammersley.F90 )
target_link_libraries( sll_random_deviate_generators
  sll_assert sll_memory sll_working_precision )

# Test
IF(BUILD_TESTING)
  add_executable( test_random_deviate_generators unit_test.F90 )
  target_link_libraries( test_random_deviate_generators
    sll_random_deviate_generators
    sll_assert
    sll_memory
    sll_working_precision )
ENDIF(BUILD_TESTING)
