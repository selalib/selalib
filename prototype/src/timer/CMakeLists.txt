set(TIMER_SRC sll_timer.F90)
add_library(sll_timer STATIC ${TIMER_SRC})
add_executable(test_timer unit_test.F90)
target_link_libraries(test_timer sll_timer sll_utilities sll_collective sll_memory sll_assert)
