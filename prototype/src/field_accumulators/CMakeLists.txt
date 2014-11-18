add_library(sll_field_accumulators STATIC sll_field_accumulators_2d.F90)
target_link_libraries(sll_field_accumulators
                      sll_memory
                      sll_assert
                      sll_working_precision)
#YG add_dependencies(sll_memory sll_assert sll_working_precision)

add_executable(test_e_field_accumulator_2d unit_test.F90)
target_link_libraries(test_e_field_accumulator_2d sll_field_accumulators)
#YG target_link_libraries(test_e_field_accumulator_2d
#YG   sll_memory 
#YG   sll_working_precision 
#YG   sll_assert
#YG   )
