if( MPI_ENABLED )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
  add_subdirectory( pic_1d )
  add_subdirectory( pic_opt2d )
  add_subdirectory( pic_basic )
  add_subdirectory( pif )
  add_subdirectory( pic_visu_parallel )
  add_subdirectory( pic_utilities )
endif( MPI_ENABLED )

ADD_SUBDIRECTORY(random_deviate_generators)
ADD_SUBDIRECTORY(pic_visu)


