#if OMPI_VERSION < 10705
  external ::       &
    mpi_bcast,      &
    mpi_allgather,  &
    mpi_allgatherv, &
    mpi_allreduce,  &
    mpi_alltoall,   &
    mpi_alltoallv,  &
    mpi_gather,     &
    mpi_gatherv,    &
    mpi_scatter,    &
    mpi_scatterv,   &
    mpi_isend,      &
    mpi_irecv,      &
    mpi_reduce,     &
    mpi_recv,       &
    mpi_send,       &
    mpi_sendrecv
#endif
