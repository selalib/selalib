ADD_EXECUTABLE( test_fmempool test_fmempool.F90 )
TARGET_LINK_LIBRARIES( test_fmempool fmempool )
ADD_TEST( NAME fmempool COMMAND test_fmempool )
SET_TESTS_PROPERTIES( fmempool PROPERTIES TIMEOUT 20 )
