ADD_EXECUTABLE( test_io test_io.F90 )
TARGET_LINK_LIBRARIES( test_io sll_file_io  )
ADD_TEST( NAME file_io COMMAND test_io )
SET_TESTS_PROPERTIES( file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_ascii_io test_ascii_io.F90 )
TARGET_LINK_LIBRARIES( test_ascii_io sll_file_io  )
ADD_TEST( NAME ascii_file_io COMMAND test_ascii_io )
SET_TESTS_PROPERTIES( ascii_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_gnuplot test_gnuplot.F90 )
TARGET_LINK_LIBRARIES( test_gnuplot sll_file_io  )
ADD_TEST( NAME gnuplot_file_io COMMAND test_gnuplot )
SET_TESTS_PROPERTIES( gnuplot_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_io_tri_mesh test_io_tri_mesh.F90 sll_generate_tri_mesh.F90 )
TARGET_LINK_LIBRARIES( test_io_tri_mesh sll_file_io )
