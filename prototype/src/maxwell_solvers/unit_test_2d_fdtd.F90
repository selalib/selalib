program test_maxwell_2d_fdtd
!------------------------------------------------------------------------
!  test 2D Maxwell solver based on finite differences on a staggered grid
!------------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use numeric_constants

use sll_maxwell
use sll_maxwell_2d_fdtd

implicit none

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

sll_int32  :: nc_eta1, nc_eta2
sll_int32  :: error

type(maxwell_fdtd)                      :: maxwell_TE
type(maxwell_fdtd)                      :: maxwell_TM

sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: bz
sll_real64, dimension(:,:), allocatable :: bz_exact

sll_real64, dimension(:,:), allocatable :: bx
sll_real64, dimension(:,:), allocatable :: by
sll_real64, dimension(:,:), allocatable :: ez
sll_real64, dimension(:,:), allocatable :: ez_exact

sll_int32                          :: i, j
sll_real64                         :: omega
sll_real64                         :: time
sll_int32                          :: istep, nstep
sll_real64                         :: err_te
sll_real64                         :: err_tm
sll_real64                         :: dt
sll_real64                         :: cfl = 0.5

sll_int32,  parameter              :: mode = 2


eta1_min = .0_f64; eta1_max = 1.0_f64
eta2_min = .0_f64; eta2_max = 1.0_f64

nc_eta1 = 127; nc_eta2 = 127

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

call initialize(maxwell_TE, 1, nc_eta1+1, 1,  &
                nc_eta2+1, delta_eta1, delta_eta2, TE_POLARIZATION)

call initialize(maxwell_TM, 1, nc_eta1+1, 1,  &
                nc_eta2+1, delta_eta1, delta_eta2, TM_POLARIZATION)

dt = cfl  / sqrt (1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))
nstep = 100

time  = 0.

omega = sqrt( (mode*sll_pi/(nc_eta1*delta_eta1))**2   &
        &    +(mode*sll_pi/(nc_eta2*delta_eta2))**2)

SLL_ALLOCATE(ex(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ey(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(bz(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(bz_exact(nc_eta1+1,nc_eta2+1), error)

SLL_ALLOCATE(bx(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(by(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ez(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ez_exact(nc_eta1+1,nc_eta2+1), error)

do istep = 1, nstep !*** Loop over time

   time = time + 0.5_f64*dt

   do j = 1, nc_eta2+1
   do i = 1, nc_eta1+1
      bz_exact(i,j) =   - cos(mode*sll_pi*(i-0.5_f64)*delta_eta1)    &
                        * cos(mode*sll_pi*(j-0.5_f64)*delta_eta2)    &
                        * cos(omega*time)
   end do
   end do

   ez_exact = - bz_exact

   if (istep == 1) then
      bz = bz_exact
      ez = ez_exact
   end if

   err_te = maxval(abs(bz_exact - bz))
   err_tm = maxval(abs(ez_exact - ez))
   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',g12.3,' sec')",advance="no") time
   write(*,"(' erreur L2 = ',2g15.5)") err_te, err_tm

   call plot_field('ez',ez, ez_exact, istep, time)

   call solve(maxwell_TE, ex, ey, bz, dt)
   call solve(maxwell_TM, bx, by, ez, dt)

   time = time + 0.5_f64*dt


end do ! next time step

print*,'PASSED'

DEALLOCATE(ex)
DEALLOCATE(ey)
DEALLOCATE(bz)
DEALLOCATE(bz_exact)

call free(maxwell_TE)

end program test_maxwell_2d_fdtd
