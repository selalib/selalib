program test_poisson_1d_pastix
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_collective
use sll_poisson_1d_pastix

implicit none

sll_real64, dimension(:), allocatable :: ex
sll_real64, dimension(:), allocatable :: ex_exact
sll_real64, dimension(:), allocatable :: rho

type (poisson_1d_pastix)            :: poisson

sll_int32   :: nc_eta1
sll_real64  :: eta1_min
sll_real64  :: eta1_max
sll_real64  :: delta_eta1
sll_real64  :: x
sll_int32   :: error
sll_int32   :: i

call sll_boot_collective()

nc_eta1 = 50

SLL_ALLOCATE(rho(nc_eta1+1),error)
SLL_ALLOCATE(ex(nc_eta1+1),error)
SLL_ALLOCATE(ex_exact(nc_eta1+1),error)

eta1_min = 0.0
eta1_max = 1.0
delta_eta1 = (eta1_max-eta1_min) / nc_eta1
do i=1,nc_eta1+1
   x = (i-1)*delta_eta1
   rho(i)      = 4*sll_pi**2*sin(2*sll_pi*x) * delta_eta1**2
   ex_exact(i) = sin(2*sll_pi*x)
   write(17,*) x, ex_exact(i)
end do

call initialize(poisson, eta1_min, eta1_max, nc_eta1, error) 

call solve(poisson, ex, rho)

do i = 1, nc_eta1+1
   write(18,*) (i-1)*delta_eta1, ex(i) 
end do
    
print*,'   error=',maxval(abs(ex-ex_exact))

call sll_halt_collective()

end program test_poisson_1d_pastix
