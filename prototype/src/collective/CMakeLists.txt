set(COLLECTIVE_SRC sll_collective.F90)
add_library(sll_collective STATIC ${COLLECTIVE_SRC})
add_executable(test_collective unit_test.F90)
target_link_libraries(test_collective sll_collective sll_timer sll_assertion  sll_memory)
