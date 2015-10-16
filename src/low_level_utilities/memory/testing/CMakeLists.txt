ADD_EXECUTABLE( test_memory test_memory.F90 )
TARGET_LINK_LIBRARIES( test_memory sll_memory )
ADD_TEST( NAME memory COMMAND test_memory )
SET_TESTS_PROPERTIES( memory PROPERTIES TIMEOUT 20 )
