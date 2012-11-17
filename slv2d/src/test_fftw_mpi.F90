program test_fftw_mpi

use fftw3

implicit none
include "mpif.h"

integer(C_INTPTR_T), parameter :: nx = 128
integer(C_INTPTR_T), parameter :: ny = 256
complex(C_DOUBLE_COMPLEX), pointer :: cdata(:,:)
integer(C_INTPTR_T) :: alloc_local
integer(C_INTPTR_T) :: alloc_local_transposed
integer(C_INTPTR_T) :: local_nx, local_i_offset
integer(C_INTPTR_T) :: local_ny, local_j_offset
integer(C_INTPTR_T) :: i, j, k
type(C_PTR) :: fw_c2c, bw_c2c, transpose_plan
type(C_PTR) :: p_cdata
integer     :: ierr, rank, task

real(8), dimension(:,:), pointer :: rho, trho

real(8)     :: tbegin, tend
integer(8)  :: ncount, t0, tn

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, task, ierr)
print *, "Hello, world, I am ", rank, " of ", task

call fftw_mpi_init()

alloc_local = fftw_mpi_local_size_2d(ny, nx, MPI_COMM_WORLD, &
              local_ny, local_j_offset)

p_cdata = fftw_alloc_complex(alloc_local)
call c_f_pointer(p_cdata, cdata, [nx,local_ny])


fw_c2c = fftw_mpi_plan_dft_2d(nx,ny,cdata,cdata, &
                              MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE)

bw_c2c = fftw_mpi_plan_dft_2d(nx,ny,cdata,cdata, &
                              MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)


call system_clock(count=t0, count_rate=ncount)
call cpu_time(tbegin)

do k = 1, 1000

   do j = 1, local_ny
      do i = 1, nx
        cdata(i, j) = k*(i+j) + local_j_offset
      end do
   end do

   call fftw_mpi_execute_dft(fw_c2c, cdata, cdata)
   call fftw_mpi_execute_dft(bw_c2c, cdata, cdata)

end do

call system_clock(count=tn,count_rate=ncount)
write(*,"(' elapsed time ', g15.5, ' s')") float(tn - t0)/real(ncount,8)
call cpu_time(tend)
write(*,"(' cpu time ', g15.5, ' s')") tend-tbegin

transpose_plan = fftw_mpi_plan_transpose(nx, ny, rho, trho, &
                 MPI_COMM_WORLD, FFTW_MEASURE)

alloc_local_transposed = fftw_mpi_local_size_2d_transposed(ny, nx, &
                         MPI_COMM_WORLD, &
                         local_ny, local_j_offset, local_nx, local_i_offset)

allocate(rho(nx,local_ny))
allocate(trho(local_nx,ny))

do k = 1, 100

   do j = 1, local_ny
      do i = 1, nx
        rho(i, j) = k*(i+j) + local_j_offset
      end do
   end do

   call fftw_execute_r2r(transpose_plan, rho, trho)
   !call fftw_execute(transpose_plan, trho, rho)

end do


call fftw_destroy_plan(fw_c2c)
call fftw_destroy_plan(bw_c2c)
call fftw_free(p_cdata)
call fftw_mpi_cleanup()
call MPI_FINALIZE(ierr)

end program test_fftw_mpi
