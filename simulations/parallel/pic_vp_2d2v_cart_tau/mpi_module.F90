module mpi_module 

use mpi
use zone

implicit none

integer, dimension(MPI_STATUS_SIZE) :: stat

integer :: nproc=1, rang=0
integer :: code, tag=1111

contains

!*********************************************************************

subroutine init_mpi( prank , psize)

integer :: prank
integer :: psize

call MPI_INIT(code)
call MPI_COMM_RANK(MPI_COMM_WORLD,prank,code)
call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code)
print*, ' Hello from mpi proc number ',prank, ' of ', psize
call MPI_BARRIER(MPI_COMM_WORLD,code)

end subroutine init_mpi

!*********************************************************************

subroutine finish_mpi()

call MPI_BARRIER(MPI_COMM_WORLD,code)
call MPI_FINALIZE(code)

end subroutine finish_mpi

!*********************************************************************

subroutine mpi_global_master()

call MPI_BCAST(nx,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
call MPI_BCAST(ny,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
call MPI_BCAST(nstep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
call MPI_BCAST(nbpart,1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
call MPI_BCAST(dimx,  1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(dimy,  1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(dx,    1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(dy,    1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(alpha, 1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(kx,    1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(ky,    1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(poids, 1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)
call MPI_BCAST(dt,    1, MPI_REAL8  , 0, MPI_COMM_WORLD, code)

end subroutine mpi_global_master

end module mpi_module
