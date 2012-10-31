program test_maxwell_2d
  !-------------------------------------------------------------------
  !  test 1D Maxwell solver based on FFT
  !-------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use numeric_constants

use sll_maxwell_2d

implicit none

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

sll_int32  :: nc_eta1, nc_eta2
sll_int32  :: error

type(maxwell_2d)                        :: maxwell_TE
sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: bz
sll_real64, dimension(:,:), allocatable :: bz_exact

sll_int32                          :: i, j
sll_real64                         :: omega
sll_real64                         :: time
sll_int32                          :: istep, nstep
sll_real64                         :: err_l2
sll_real64                         :: dt
sll_real64                         :: cfl = 0.5

sll_int32,  parameter              :: mode = 2

character(len=4)                   :: counter
sll_real64, dimension(:,:), allocatable :: x
sll_real64, dimension(:,:), allocatable :: y

sll_real64, dimension(:,:), allocatable :: hx
sll_real64, dimension(:,:), allocatable :: hy
sll_real64, dimension(:,:), allocatable :: ez

!sll_real64, dimension(:,:), allocatable :: ex
!sll_real64, dimension(:,:), allocatable :: ey
!sll_real64, dimension(:,:), allocatable :: hz
sll_real64 :: xmin, xmax, ymin, ymax
sll_real64 :: omega, dt, time, dx, dy
sll_int32  :: nstep, istep, iplot, error
sll_int32  :: nx, ny, i, j, ik, jk

sll_real64, dimension(:,:), allocatable :: tmp

type(maxwell_pstd) :: this

nx = 128
ny = 128

SLL_ALLOCATE(hx(nx,ny), error)
SLL_ALLOCATE(hy(nx,ny), error)
SLL_ALLOCATE(ez(nx,ny), error)

!SLL_ALLOCATE(ex(nx,ny), error)
!SLL_ALLOCATE(ey(nx,ny), error)
!SLL_ALLOCATE(hz(nx,ny), error)

SLL_ALLOCATE(tmp(nx,ny), error)

xmin = 0.
xmax = sll_pi
ymin = 0.
ymax = sll_pi

call init_maxwell_solver(this, xmin, xmax, nx, &
                               ymin, ymax, ny, &
                               ez, error )                                

SLL_ALLOCATE(x(nx,ny),error)
SLL_ALLOCATE(y(nx,ny),error)
dx =  (xmax - xmin) / (nx-1)
dy =  (ymax - ymin) / (ny-1)

do j = 1, ny
   do i = 1, nx
      x(i,j) = xmin + i*dx
      y(i,j) = ymin + j*dy
   end do
end do

dt = 0.5_f64  / sqrt (1./(dx*dx)+1./(dy*dy)) 
nstep = floor(1./dt)
write(*,*)
write(*,*) " dx = ", dx
write(*,*) " dy = ", dy
write(*,*) " dt = ", dt
write(*,*)
write(*,*) " Nombre d'iteration nstep = ", nstep

time  = 0.
iplot = 0

omega = sqrt(2._f64)

!Polarisation TE
!ex =  cos(x)*sin(y)*sin(omega*time)/omega
!ey = -sin(x)*cos(y)*sin(omega*time)/omega
!hz = -cos(x)*cos(y)*cos(omega*time)
!
!Polarisation TM
!ez =  cos(x)*cos(y)*cos(omega*time)
!hx =  cos(x)*sin(y)*sin(omega*time)/omega
!hy = -sin(x)*cos(y)*sin(omega*time)/omega

time = 0._f64
hx = 0.0
hy = 0.0
ez = cos(x)*cos(y)

do istep = 1, nstep ! Loop over time ************************

   call solve_faraday(this, hx, hy, ez)
   
   call solve_ampere_maxwell(this, hx, hy, ez)

   time = time + dt
   write(*,*) istep , maxval(abs(ez - cos(x)*cos(y)*cos(omega*time)))

end do ! Next time step *************************************

eta1_min = .0_f64; eta1_max = 1.0_f64
eta2_min = .0_f64; eta2_max = 1.0_f64

nc_eta1 = 127; nc_eta2 = 127

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

call new(maxwell_TE, 1, nc_eta1+1, 1, nc_eta2+1, delta_eta1, delta_eta2)

dt = cfl  / sqrt (1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))
nstep = 100

time  = 0.

omega = sqrt( (mode*sll_pi/(nc_eta1*delta_eta1))**2   &
        &    +(mode*sll_pi/(nc_eta2*delta_eta2))**2)

SLL_ALLOCATE(ex(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ey(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(bz(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(bz_exact(nc_eta1+1,nc_eta2+1), error)

do istep = 1, nstep !*** Loop over time

   time = time + 0.5_f64*dt

   do i=1,nc_eta1+1
   do j=1,nc_eta2+1
      bz_exact(i,j) =   - cos(mode*sll_pi*(i-0.5_f64)/nc_eta1)    &
                        * cos(mode*sll_pi*(j-0.5_f64)/nc_eta2)    &
                        * cos(omega*time)
   end do  
   end do  

   if (istep == 1) bz = bz_exact

   call solve(maxwell_TE, ex, ey, bz, dt)

   time = time + 0.5_f64*dt

   err_l2 = maxval(abs(bz_exact - bz))
   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',g12.3,' sec')",advance="no") time
   write(*,"(' erreur L2 = ',g10.5)") sqrt(err_l2)

   call int2string(istep, counter)

end do ! next time step

print*,'PASSED'

DEALLOCATE(ex)
DEALLOCATE(ey)
DEALLOCATE(bz)
DEALLOCATE(bz_exact)

end program test_maxwell_2d
