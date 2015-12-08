#if OMPI_V_MINOR < 10
  external :: &
#if OMPI_V_MINOR < 8
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
#endif
