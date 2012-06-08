add_library(sll_assertion STATIC sll_assert.F90)
add_executable(test_assert unit_test.F90)
target_link_libraries(test_assert sll_assertion)
