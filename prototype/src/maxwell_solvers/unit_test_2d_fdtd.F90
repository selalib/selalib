program test_maxwell_2d_pstd
  !-------------------------------------------------------------------
  !  test 1D Maxwell solver based on FFT
  !-------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use numeric_constants

use sll_maxwell_2d_pstd

implicit none

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

sll_int32  :: nc_eta1, nc_eta2
sll_int32  :: error

type(maxwell_pstd)                      :: maxwell_TM
sll_real64, dimension(:,:), allocatable :: hx
sll_real64, dimension(:,:), allocatable :: hy
sll_real64, dimension(:,:), allocatable :: ez

sll_int32                          :: i, j
sll_real64                         :: omega
sll_real64                         :: time
sll_int32                          :: istep, nstep
sll_real64                         :: err_l2
sll_real64                         :: dt
sll_real64                         :: cfl = 0.5

sll_int32,  parameter              :: mode = 2

character(len=4)                   :: counter

eta1_min = .0_f64; eta1_max = 1.0_f64
eta2_min = .0_f64; eta2_max = 1.0_f64

nc_eta1 = 127; nc_eta2 = 127

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

call initialize(maxwell_TM, eta1_min, eta1_max, nc_eta1+1, &
                            eta2_min, eta2_max, nc_eta2+1, ez, error)

dt = cfl  / sqrt (1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))
nstep = 100

time  = 0.

omega = sqrt( (mode*sll_pi/(nc_eta1*delta_eta1))**2   &
        &    +(mode*sll_pi/(nc_eta2*delta_eta2))**2)

SLL_ALLOCATE(hx(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(hy(nc_eta1+1,nc_eta2+1), error)
SLL_ALLOCATE(ez(nc_eta1+1,nc_eta2+1), error)

do istep = 1, nstep !*** Loop over time

   time = time + 0.5_f64*dt

      ez_exact(i,j) =   - cos(mode*sll_pi*(i-0.5_f64)/nc_eta1)    &
                        * cos(mode*sll_pi*(j-0.5_f64)/nc_eta2)    &
                        * cos(omega*time)

   if (istep == 1) bz = bz_exact

   call solve(maxwell_TM, ex, ey, bz, dt)

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

call free(maxwell_TM)

end program test_maxwell_2d_pstd
