program hello
use mpi
use petsc_module
implicit none
integer :: prank, psize, code
integer :: nthreads
integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

call MPI_INIT(code)                             !---- Initialisation MPI
call MPI_Comm_rank(MPI_COMM_WORLD,prank,code)
call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code)

call hello_petsc()
  
!Fork a team of threads giving them their own copies of variables
!$OMP PARALLEL PRIVATE(NTHREADS, TID)
!Obtain thread number
!$OMP CRITICAL
PRINT *, psize,prank,OMP_GET_NUM_THREADS(),OMP_GET_THREAD_NUM()
!$OMP END CRITICAL
!$OMP PARALLEL

call MPI_FINALIZE(code)

end program hello
