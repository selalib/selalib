add_subdirectory(file_io)
add_subdirectory(xdmf_io) # YG [work in progress]

if(MPI_ENABLED AND HDF5_IS_PARALLEL)

  add_subdirectory(file_io_parallel)
  add_subdirectory(xdmf_io_parallel) # YG [work in progress]

endif(MPI_ENABLED AND HDF5_IS_PARALLEL)
