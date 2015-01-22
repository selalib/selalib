!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
!-------------------------------------------------------------------
!  test 1D Ampere solver based on FFT
!-------------------------------------------------------------------
program test_ampere_1d_pstd
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_module_ampere_1d_pstd

implicit none

sll_real64 :: x_max, x_min, delta_x
sll_int32  :: nc_x
sll_int32  :: error

type(sll_ampere_1d_pstd) :: ampere
sll_int32                :: i
sll_real64               :: w
sll_real64               :: time
sll_int32                :: istep, nstep
sll_real64               :: dt
sll_real64               :: cfl = 0.5

sll_int32,  parameter    :: mode = 2

sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: ex
sll_real64, dimension(:), allocatable :: ex_exact
sll_real64, dimension(:), allocatable :: jx
sll_real64, dimension(:), allocatable :: rho

sll_real64 :: tstart, tend

call cpu_time(tstart)

x_min = .0_f64; x_max = 1.0_f64

nc_x = 127

delta_x = (x_max-x_min)/nc_x

SLL_ALLOCATE(x(nc_x+1),error)

do i = 1, nc_x+1
  x(i) = x_min + (i-1)*delta_x
end do

dt = cfl  / sqrt(1./(delta_x*delta_x))
nstep = 100

time  = 0.

w = 2 * sll_pi

SLL_ALLOCATE(ex(nc_x+1), error)
SLL_ALLOCATE(jx(nc_x+1), error)
SLL_ALLOCATE(ex_exact(nc_x+1), error)
SLL_ALLOCATE(rho(nc_x+1), error)

call sll_create(ampere, x_min, x_max, nc_x)

do istep = 1, nstep !*** Loop over time

  ex  = w*sin(time*w)*sin(sll_pi*x)
  rho = -w*w*cos(time*w)*sin(sll_pi*x)
  jx  = sll_pi*w*cos(sll_pi*x)*sin(time*w)
  ex_exact = ex
  time = time + 0.5_f64*dt

  write(*,"(10x,' istep = ',I6)",advance="no") istep
  write(*,"(' time = ',g12.3,' sec')",advance="no") time
  write(*,"(' erreurs TM,TE = ',g15.5)") maxval(abs(ex - ex_exact))

  call sll_solve(ampere, ex, dt, jx)

  time = time + 0.5_f64*dt

end do ! next time step

call cpu_time(tend)
print"('CPU time : ',g15.3)", tend-tstart
print*,'PASSED'

DEALLOCATE(ex)
DEALLOCATE(jx)
DEALLOCATE(rho)
DEALLOCATE(ex_exact)

end program test_ampere_1d_pstd
