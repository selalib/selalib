ADD_EXECUTABLE( test_utilities test_utilities.F90 )
TARGET_LINK_LIBRARIES( test_utilities sll_utilities )
ADD_TEST( NAME utilities COMMAND test_utilities )

ADD_EXECUTABLE( test_fornberg test_fornberg.F90 )
TARGET_LINK_LIBRARIES( test_fornberg sll_utilities )
ADD_TEST( NAME fornberg COMMAND test_fornberg )
