program test_poisson_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
  use sll_cubic_splines
  use sll_poisson_1D_periodic

  implicit none

  sll_real64, dimension(:), allocatable :: ex
  sll_real64, dimension(:), allocatable :: ex_exact
  sll_real64, dimension(:), allocatable :: rho

  type (poisson_1d_periodic)            :: poisson

  sll_int32   :: nc_eta1
  sll_real64  :: eta1_min
  sll_real64  :: eta1_max
  sll_real64  :: delta_eta1
  sll_real64  :: x
  sll_int32   :: error
  sll_int32   :: mode
  sll_int32   :: i

  nc_eta1 = 128

  SLL_ALLOCATE(rho(nc_eta1+1),error)
  SLL_ALLOCATE(ex(nc_eta1+1),error)
  SLL_ALLOCATE(ex_exact(nc_eta1+1),error)

  eta1_min = 0.0
  eta1_max = 2*sll_pi
  delta_eta1 = (eta1_max-eta1_min) / nc_eta1
  mode = 4
  do i=1,nc_eta1+1
     x = (i-1)*delta_eta1
     rho(i)      =  mode**2*sin(mode*x)
     ex_exact(i) = -mode*cos(mode*x)
  end do

  call initialize(poisson, eta1_min, eta1_max, nc_eta1, error) 

  call solve(poisson, ex, rho)

  print*,'mode=',mode,'   error=',maxval(abs(ex-ex_exact))

  if (error < 1.e-14) then
     print*, 'PASSED'
  else 
     print*, 'FAILED'
     stop
  end if


end program test_poisson_1d
