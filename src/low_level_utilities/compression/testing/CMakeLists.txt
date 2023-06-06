ADD_EXECUTABLE( test_compression test_compression.F90 )
TARGET_LINK_LIBRARIES( test_compression sll_compression )
ADD_TEST( NAME compression COMMAND test_compression )
#SET_TESTS_PROPERTIES( compression PROPERTIES TIMEOUT 20 )
