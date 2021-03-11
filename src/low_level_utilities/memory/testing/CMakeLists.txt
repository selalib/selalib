add_executable(test_memory test_memory.F90)
target_link_libraries(test_memory sll_memory)
add_test(NAME memory COMMAND test_memory)
set_tests_properties(memory PROPERTIES TIMEOUT 20)
