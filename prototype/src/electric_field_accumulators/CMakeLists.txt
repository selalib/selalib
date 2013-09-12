add_library(sll_electric_field_accumulators STATIC sll_electric_field_2d.F90)
add_dependencies(sll_memory sll_assertion sll_working_precision)

add_executable(test_e_field_accumulator_2d unit_test.F90)
target_link_libraries(test_e_field_accumulator_2d
  sll_memory 
  sll_working_precision 
  sll_assertion
  )
