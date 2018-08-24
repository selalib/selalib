ADD_SUBDIRECTORY( file_io )
ADD_SUBDIRECTORY( xdmf_io ) # YG [work in progress]

IF(MPI_ENABLED AND HDF5_IS_PARALLEL)
   
  ADD_SUBDIRECTORY( file_io_parallel )
  ADD_SUBDIRECTORY( xdmf_io_parallel ) # YG [work in progress]

ENDIF(MPI_ENABLED AND HDF5_IS_PARALLEL)
