# setting libraries for GMP and Boost Test Framework
# list of sources used in the kernel part of the lib
set(MEMORY_SRC sll_memory.F90)
# bigint test
add_library(sll_memory STATIC ${MEMORY_SRC})
add_executable(test_memory unit_test.F90)
target_link_libraries(test_memory sll_memory)
