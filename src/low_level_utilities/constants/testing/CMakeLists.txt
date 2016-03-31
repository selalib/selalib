ADD_EXECUTABLE( test_constants test_constants.F90 )
ADD_TEST( NAME constants COMMAND test_constants )
TARGET_LINK_LIBRARIES( test_constants sll_constants )
