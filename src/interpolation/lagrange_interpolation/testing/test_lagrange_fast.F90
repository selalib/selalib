program test_lagrange_fast
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
#include "sll_working_precision.h"
  
  use sll_m_constants, only : sll_pi
  
  use sll_m_lagrange_fast, only : &
    sll_s_interpolate_array_disp_lagrange_fixed_no_bc
  

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32 :: i, s, num_points
  sll_real64 :: diff, alpha, xmin, xmax, l
  sll_real64, dimension(:), allocatable :: xi, fi, xp, fp

  s = 5  ! stencil
  num_points = 100
  alpha = 0.2_f64
  ! known function values
  allocate(xi(1:num_points))
  allocate(fi(1:num_points))
  ! x values with offset
  allocate(xp(1:num_points))
  ! interpolated function values
  allocate(fp(1:num_points))

  ! data initialization
  xmin = 0.0_f64
  xmax = num_points-1.0_f64
  l = xmax - xmin
  do i = 1, num_points
    xi(i) = real(i - 1, f64)
    fi(i) = f(xi(i), num_points)
    xp(i) = xi(i) + alpha
  end do
  diff = 0.0_f64
  fp(:) = 0.0_f64

  call sll_s_interpolate_array_disp_lagrange_fixed_no_bc(fi, fp, alpha, s)

  do i=1,num_points
    diff=max(diff, abs(f(xp(i),num_points) - fp(i)))
  end do

  if (diff < 1.e-4) then
    print *, ""
    print *, "Fast Lagrange interpolation unit test: PASSED"
    print *, "error =", diff
  else
    print *, ""
    print *, "Fast Lagrange interpolation unit test: FAILED"
    print *, "error =", diff
  end if

  ! TODO: implement tests for periodic and halo routines

  deallocate (xi,fi,xp,fp)
  stop

contains

  function f(x, num_points)
    sll_int32, intent(in) :: num_points
    sll_real64, intent(in) :: x
    sll_real64 :: f
    f = cos(2*sll_pi*x/(num_points-1.0))
  end function

end program test_lagrange_fast
