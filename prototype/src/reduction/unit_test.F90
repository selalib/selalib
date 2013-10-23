program test_reduction
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_reduction_module
implicit none

  sll_real64, dimension(:,:,:,:), allocatable   :: data_4d
  sll_real64, dimension(:,:,:), allocatable   :: data_3d
  sll_int32 :: ierr
  sll_real64 :: err
  sll_int32 :: Npts(4) = (/33,33,33,65/)
  sll_real64 :: delta
  
  SLL_ALLOCATE(data_4d(Npts(1),Npts(2),Npts(3),Npts(4)),ierr)
  SLL_ALLOCATE(data_3d(Npts(1),Npts(2),Npts(3)),ierr)
  
  
  delta = 1._f64/real(Npts(4)-1,f64)
  
  data_4d = 1._f64
  
  call compute_reduction_4d_to_3d_direction4(&
    data_4d, &
    data_3d, &
    Npts(1), &
    Npts(2), &
    Npts(3), &
    Npts(4), &
    delta)

  err=maxval(abs(data_3d-1._f64))
  
  
  print *,'#err=',err

  if(err==0._f64)then
    print *,'#PASSED'
  endif

end program test_reduction