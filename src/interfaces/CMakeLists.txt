add_subdirectory( sparse_matrix_manager )
add_subdirectory( fft )

if( MPI_ENABLED )
  add_subdirectory( mpi )
endif( MPI_ENABLED )
