add_library(sll_timer STATIC sll_timer.F90)
add_executable(test_timer unit_test.F90)
target_link_libraries(test_timer sll_timer)
