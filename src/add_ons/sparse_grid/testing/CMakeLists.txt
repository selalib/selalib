ADD_EXECUTABLE( test_sparse_grid_2d test_sparse_grid_2d.F90 )
TARGET_LINK_LIBRARIES (test_sparse_grid_2d sll_sparse_grid sll_interpolators sll_periodic_interpolation)
ADD_TEST( NAME sparse_grid_2d  COMMAND test_sparse_grid_2d )
SET_TESTS_PROPERTIES( sparse_grid_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_sparse_grid_3d test_sparse_grid_3d.F90 )
TARGET_LINK_LIBRARIES (test_sparse_grid_3d sll_sparse_grid sll_interpolators sll_periodic_interpolation)
ADD_TEST( NAME sparse_grid_3d  COMMAND test_sparse_grid_3d )
SET_TESTS_PROPERTIES( sparse_grid_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_sparse_grid_4d test_sparse_grid_4d.F90 )
TARGET_LINK_LIBRARIES (test_sparse_grid_4d sll_sparse_grid sll_interpolators sll_periodic_interpolation)
ADD_TEST( NAME sparse_grid_4d  COMMAND test_sparse_grid_4d )
SET_TESTS_PROPERTIES( sparse_grid_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_2d_sparse_grid_fft test_poisson_2d_sparse_grid_fft.F90 )
TARGET_LINK_LIBRARIES( test_poisson_2d_sparse_grid_fft 
      ${POISSON_LIBS} sll_sparse_grid ${FFTW_LIBRARIES} )
ADD_TEST( NAME poisson_2d_sparse_grid_fft  COMMAND test_poisson_2d_sparse_grid_fft )
SET_TESTS_PROPERTIES( poisson_2d_sparse_grid_fft PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_poisson_3d_sparse_grid_fft test_poisson_3d_sparse_grid_fft.F90 )
TARGET_LINK_LIBRARIES( test_poisson_3d_sparse_grid_fft 
  ${POISSON_LIBS} sll_sparse_grid ${FFTW_LIBRARIES} )
ADD_TEST( NAME poisson_3d_sparse_grid_fft  COMMAND test_poisson_3d_sparse_grid_fft )
SET_TESTS_PROPERTIES( poisson_3d_sparse_grid_fft PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )