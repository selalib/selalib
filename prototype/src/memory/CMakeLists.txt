set(MEMORY_SRC sll_memory.F90)
add_library(sll_memory STATIC ${MEMORY_SRC})
add_executable(test_memory unit_test.F90)
target_link_libraries(test_memory sll_memory)
