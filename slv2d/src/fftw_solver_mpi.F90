program fftw_solver_mpi 
use fftw3

#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_constants

implicit none
include "mpif.h"

integer(C_INTPTR_T), parameter :: nx = 128
integer(C_INTPTR_T), parameter :: ny = 256
integer(C_INTPTR_T) :: alloc_local_rhot
integer(C_INTPTR_T) :: local_nx, local_i_offset
integer(C_INTPTR_T) :: local_ny, local_j_offset
integer(C_INTPTR_T) :: i, j, k
type(C_PTR) :: fw_r2c, bw_c2r
integer     :: ierr, rank, task

sll_int32   :: error
sll_real64  :: xmax, xmin, ymax, ymin
sll_real64  :: x, y, delta_x, delta_y

sll_real64, dimension(:,:), pointer :: phi, phi_exact
sll_int32                           :: mode
sll_real64                          :: kx0, kx
sll_real64                          :: ky0, ky
sll_int32                           :: ik, jk

sll_real64, dimension(:,:), allocatable :: kmod
sll_real64, dimension(:,:), pointer :: rho
complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: rhot
type(C_PTR) :: p_rhot

sll_real64 :: tbegin, tend
sll_int64  :: ncount, t0, tn

sll_int32, parameter :: nthreads = 4

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, task, ierr)
print *, "Hello, world, I am ", rank, " of ", task

call fftw_mpi_init()

xmin = .0_f64; xmax = 2.0_f64*sll_pi
ymin = .0_f64; ymax = 2.0_f64*sll_pi

delta_x = (xmax-xmin) / nx
delta_y = (ymax-ymin) / ny

kx0=2._f64*sll_pi/(xmax-xmin)
ky0=2._f64*sll_pi/(ymax-ymin)

SLL_ALLOCATE(kmod(nx/2+1,ny), error)

do ik=1,nx/2+1
   kx  = (ik-1)*kx0
   do jk = 1, ny/2
      ky  = (jk-1)*ky0
      kmod(ik,jk) = kx*kx+ky*ky
   end do
   do jk = ny/2+1,ny     
      ky= (jk-1-ny)*ky0
      kmod(ik,jk) = kx*kx+ky*ky
   end do
end do
kmod(1,1) = 1.0_f64

alloc_local_rhot = fftw_mpi_local_size_2d(ny, nx, MPI_COMM_WORLD, &
              local_ny, local_j_offset)

SLL_ALLOCATE(rho(nx,local_ny), error)
SLL_ALLOCATE(phi(nx,local_ny), error)
SLL_ALLOCATE(phi_exact(nx,local_ny), error)

mode = 2
do i = 1, nx
   do j = 1, local_ny
      x = xmin+(i-1)*delta_x
      y = ymin+(j-1)*delta_y
      phi(i,j) = mode * sin(mode*x) * cos(mode*y)
      rho(i,j) = 2_f64 * mode**3 * sin(mode*x)*cos(mode*y)
   end do
end do

phi_exact = phi
phi = 0.0

p_rhot = fftw_alloc_complex(alloc_local_rhot)
call c_f_pointer(p_rhot, rhot, [nx,local_ny])

fw_r2c = fftw_mpi_plan_dft_r2c_2d(nx,ny,rho,rhot, &
                              MPI_COMM_WORLD,FFTW_MEASURE)

bw_c2r = fftw_mpi_plan_dft_c2r_2d(nx,ny,rhot,rho, &
                              MPI_COMM_WORLD,FFTW_MEASURE)

call system_clock(count=t0, count_rate=ncount)
call cpu_time(tbegin)

do k = 1, 10

   call fftw_mpi_execute_dft_r2c(fw_r2c, rho, rhot)
   rhot = rhot / kmod
   call fftw_mpi_execute_dft_c2r(bw_c2r, rhot, phi)

   write(*,*) " E = ", maxval(phi_exact-phi)

end do

call system_clock(count=tn,count_rate=ncount)
write(*,"(' elapsed time ', g15.5, ' s')") float(tn - t0)/real(ncount,f64)
call cpu_time(tend)
write(*,"(' cpu time ', g15.5, ' s')") tend-tbegin

call fftw_free(p_rhot)
call fftw_mpi_cleanup()
call MPI_FINALIZE(ierr)

end program fftw_solver_mpi
