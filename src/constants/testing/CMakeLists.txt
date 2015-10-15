ADD_EXECUTABLE( test_constants test_constants.F90 )
TARGET_LINK_LIBRARIES(test_constants sll_constants)
ADD_TEST( NAME constants COMMAND test_constants )
