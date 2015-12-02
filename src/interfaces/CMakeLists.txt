add_subdirectory( sparse_matrix_manager )
add_subdirectory( fft )

if( MPI_ENABLED )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_Compiler} )
  add_subdirectory( mpi )
endif( MPI_ENABLED )
