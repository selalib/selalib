set( CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS_RELEASE} )

add_library( nufft STATIC
  dirft1d.f dirft2d.f next235.f nufft1df90.f nufft2df90.f dfftpack.f )

install( TARGETS nufft ARCHIVE DESTINATION lib )

# Tests
if( BUILD_TESTING )
  # 1D test
  add_executable( test_nufft1d nufft1d_demof90.f )
  target_link_libraries( test_nufft1d nufft )
  # 2D test
  add_executable( test_nufft2d nufft2d_demof90.f )
  target_link_libraries( test_nufft2d nufft )
endif()
