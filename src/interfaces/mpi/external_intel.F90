  external ::        &
    mpi_allgather,   &
    mpi_allgatherv,  &
    mpi_allreduce,   &
    mpi_alltoall,    &
    mpi_alltoallv,   &
    mpi_bcast,       &
    mpi_cart_coords, &
    mpi_cart_create, &
    mpi_cart_get,    &
    mpi_dims_create, &
#if INTEL_MPI_VERSION < 50000
    mpi_finalize,    &
#endif
    mpi_gather,      &
    mpi_gatherv,     &
    mpi_iallreduce,  &
    mpi_irecv,       &
    mpi_isend,       &
    mpi_reduce,      &
    mpi_recv,        &
    mpi_scatter,     &
    mpi_scatterv,    &
    mpi_send,        &
    mpi_sendrecv
