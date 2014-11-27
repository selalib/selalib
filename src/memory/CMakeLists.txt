SET(MEMORY_SRC sll_memory.F90)
ADD_LIBRARY(sll_memory STATIC ${MEMORY_SRC})
ADD_EXECUTABLE(test_memory unit_test.F90)
TARGET_LINK_LIBRARIES(test_memory sll_memory)

# Ctest
ADD_TEST( NAME memory COMMAND test_memory )
SET_TESTS_PROPERTIES( memory PROPERTIES TIMEOUT 20 )
