add_executable( test_random_deviate_generators test_random_deviate_generators.F90 )
target_link_libraries( test_random_deviate_generators
  sll_random_deviate_generators
  sll_assert
  sll_memory
  sll_working_precision )

