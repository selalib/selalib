add_library(sll_field_accumulators STATIC sll_field_accumulators_2d.F90)
target_link_libraries(sll_field_accumulators
                      sll_memory
                      sll_assert
                      sll_working_precision)

# Ctest
IF(BUILD_TESTING)
  add_executable( test_e_field_accumulator_2d unit_test.F90 )
  target_link_libraries( test_e_field_accumulator_2d sll_field_accumulators )
  ADD_TEST( NAME electric_field_accumulators COMMAND test_e_field_accumulator_2d )
  SET_TESTS_PROPERTIES( electric_field_accumulators PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)
