program test_maxwell_2d_pstd
!-------------------------------------------------------------------
!  test 2D Maxwell solver based on FFT
!-------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_maxwell_solvers.h"
use sll_constants

use sll_maxwell_2d_pstd

implicit none

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

sll_int32  :: nc_eta1, nc_eta2
sll_int32  :: error

type(maxwell_2d_pstd)              :: maxwell_TE
type(maxwell_2d_pstd)              :: maxwell_TM
sll_int32                          :: i, j
sll_real64                         :: omega
sll_real64                         :: time
sll_int32                          :: istep, nstep
sll_real64                         :: err_tm, err_te
sll_real64                         :: dt
sll_real64                         :: cfl = 0.5

sll_int32,  parameter              :: mode = 2

sll_real64, dimension(:,:), allocatable :: eta1
sll_real64, dimension(:,:), allocatable :: eta2

sll_real64, dimension(:,:), allocatable :: hx
sll_real64, dimension(:,:), allocatable :: hy
sll_real64, dimension(:,:), allocatable :: ez, ez_exact

sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: hz, hz_exact
sll_real64 :: tstart, tend

call cpu_time(tstart)
!Polarisation TE
!ex =  cos(x)*sin(y)*sin(omega*time)/omega
!ey = -sin(x)*cos(y)*sin(omega*time)/omega
!hz = -cos(x)*cos(y)*cos(omega*time)
!
!Polarisation TM
!ez =  cos(x)*cos(y)*cos(omega*time)
!hx =  cos(x)*sin(y)*sin(omega*time)/omega
!hy = -sin(x)*cos(y)*sin(omega*time)/omega

eta1_min = .0_f64; eta1_max = 1.0_f64
eta2_min = .0_f64; eta2_max = 1.0_f64

nc_eta1 = 127; nc_eta2 = 127

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

SLL_ALLOCATE(eta1(nc_eta1+1,nc_eta2+1),error)
SLL_ALLOCATE(eta2(nc_eta1+1,nc_eta2+1),error)

do j = 1, nc_eta2+1
   do i = 1, nc_eta1+1
      eta1(i,j) = eta1_min + (i-1)*delta_eta1
      eta2(i,j) = eta2_min + (j-1)*delta_eta2
   end do
end do

dt = cfl  / sqrt (1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))
nstep = 100

time  = 0.

omega = sqrt( (mode*sll_pi/(nc_eta1*delta_eta1))**2   &
        &    +(mode*sll_pi/(nc_eta2*delta_eta2))**2)

SLL_ALLOCATE(hx(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(hy(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ez(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ez_exact(nc_eta2+1,nc_eta2+1), error)

call initialize(maxwell_TM, eta1_min, eta1_max, nc_eta1, &
         eta2_min, eta2_max, nc_eta2, TM_POLARIZATION)

SLL_ALLOCATE(ex(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ey(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(hz(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(hz_exact(nc_eta2+1,nc_eta2+1), error)

call initialize(maxwell_TE, eta1_min, eta1_max, nc_eta1, &
         eta2_min, eta2_max, nc_eta2, TE_POLARIZATION)


do istep = 1, nstep !*** Loop over time

   time = time + 0.5_f64*dt

   ez_exact =   cos(mode*sll_pi*eta1)*cos(mode*sll_pi*eta2)*cos(omega*time)
   hz_exact = - cos(mode*sll_pi*eta1)*cos(mode*sll_pi*eta2)*cos(omega*time)
  
   if (istep == 1) then
      ez = ez_exact
      hz = hz_exact
   end if

   !call plot_fields('ez',ez, ez_exact, istep, time)
   !call plot_fields('hz',hz, hz_exact, istep, time)

   err_tm = maxval(abs(ez - ez_exact))
   err_te = maxval(abs(hz - hz_exact))

   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',g12.3,' sec')",advance="no") time
   write(*,"(' erreurs TM,TE = ',2g15.5)") err_tm, err_te

   call solve(maxwell_TM, hx, hy, ez, dt)

   call solve(maxwell_TE, ex, ey, hz, dt)

   time = time + 0.5_f64*dt

end do ! next time step

call cpu_time(tend)
print"('CPU time : ',g15.3)", tend-tstart
print*,'PASSED'

DEALLOCATE(hx)
DEALLOCATE(hy)
DEALLOCATE(ez)

end program test_maxwell_2d_pstd
