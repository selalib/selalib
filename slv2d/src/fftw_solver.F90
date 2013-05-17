program fftw_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
#include "sll_constants.h"
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

type(C_PTR) :: fw, bw
sll_int32   :: nx, ny
sll_int32   :: error
sll_int32   :: nc_eta1, nc_eta2
sll_real64  :: eta1_max, eta1_min, eta2_max, eta2_min
sll_real64  :: delta_eta1, delta_eta2
sll_real64, dimension(:,:), pointer :: phi, phi_exact, rho
sll_real64                          :: x1, x2
sll_int32                           :: mode
sll_int32                           :: i, j
sll_real64                          :: dx,dy
sll_real64                          :: kx0, kx
sll_real64                          :: ky0, ky
sll_int32                           :: ik, jk
sll_int32                           :: k

sll_real64, dimension(:,:), allocatable :: kmod

complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: rhot
integer(C_SIZE_T) :: nxh1
type(C_PTR) :: p_rhot

sll_real64 :: tbegin, tend
sll_int64  :: ncount, t0, tn

sll_int32, parameter :: nthreads = 4

eta1_min = .0_f64; eta1_max = 2.0_f64*sll_pi
eta2_min = .0_f64; eta2_max = 2.0_f64*sll_pi

nc_eta1 = 1024; nc_eta2 = 512

delta_eta1 = (eta1_max-eta1_min) / nc_eta1
delta_eta2 = (eta2_max-eta2_min) / nc_eta2

SLL_ALLOCATE(rho(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(phi(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(phi_exact(nc_eta1+1,nc_eta2+1), error)

mode = 2
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      x1 = eta1_min+(i-1)*delta_eta1
      x2 = eta2_min+(j-1)*delta_eta2
      phi(i,j) = mode * sin(mode*x1) * cos(mode*x2)
      rho(i,j) = 2_f64 * mode**3 * sin(mode*x1)*cos(mode*x2)
   end do
end do

phi_exact = phi
phi = 0.0

nx = nc_eta1
ny = nc_eta2

dx = delta_eta1
dy = delta_eta2

nxh1 = int((nx/2+1)*ny,C_SIZE_T)
p_rhot = fftw_alloc_complex(nxh1)
call c_f_pointer(p_rhot, rhot, [nx/2+1,ny])

#ifdef FFTW_THREADS
call dfftw_init_threads(error)
if (error == 0) stop 'FFTW CAN''T USE THREADS'
call dfftw_plan_with_nthreads(nthreads)
#endif 

fw = fftw_plan_dft_r2c_2d(ny, nx, rho(1:nx,1:ny), rhot, FFTW_PATIENT)
bw = fftw_plan_dft_c2r_2d(ny, nx, rhot, phi(1:nx,1:ny), FFTW_PATIENT)

kx0=2._f64*sll_pi/(eta1_max-eta1_min)
ky0=2._f64*sll_pi/(eta2_max-eta2_min)

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

call system_clock(count=t0, count_rate=ncount)
call cpu_time(tbegin)


do k = 1, 100

   call fftw_execute_dft_r2c(fw, rho(1:nx,1:ny), rhot)

   rhot = rhot / kmod

   call fftw_execute_dft_c2r(bw,rhot,phi(1:nx,1:ny))

   phi = phi / (nx*ny)

   phi(nx+1,:) = phi(1,:)
   phi(:,ny+1) = phi(:,1)

   write(*,*) " E = ", maxval(phi_exact-phi)

end do
call system_clock(count=tn,count_rate=ncount)
write(*,"(' elapsed time ', g15.5, ' s')") float(tn - t0)/ncount
call cpu_time(tend)
write(*,"(' cpu time ', g15.5, ' s')") tend-tbegin


#ifdef FFTW_THREADS
call dfftw_cleanup_threads(error)
#endif
call fftw_free(p_rhot)
call fftw_destroy_plan(fw)
call fftw_destroy_plan(bw)

end program fftw_solver
