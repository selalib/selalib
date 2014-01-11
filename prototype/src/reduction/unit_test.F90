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
  sll_real64, external :: my_compute_integral_trapezoid_1d
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
  
  call compute_reduction_4d_to_3d_direction4(&
    data_4d, &
    data_3d, &
    Npts(1), &
    Npts(2), &
    Npts(3), &
    Npts(4), &
    delta, &
    my_compute_integral_trapezoid_1d)
  print *,'#err1=',err

  err=err+maxval(abs(data_3d-1._f64))


  
  print *,'#err2=',err

  if(err==0._f64)then
    print *,'#PASSED'
  endif

end program test_reduction



  function my_compute_integral_trapezoid_1d(data, Npts, delta, func_params) &
    result(res)

    real(8), dimension(:), intent(in)    :: data
    integer, intent(in) :: Npts
    real(8),intent(in) :: delta
    real(8), dimension(:), intent(in) ,optional :: func_params
    real(8) :: res
    integer :: i

    if(present(func_params))then
      if(size(func_params)>100)then
        print *,'#size of func_params is >100'
      endif
    endif

    
    res = 0.5*(data(1)+data(Npts))
    do i=2,Npts-1
      res = res + data(i)
    enddo
    res = res*delta
  end function my_compute_integral_trapezoid_1d
