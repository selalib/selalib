program test_lagrange_interpolation_1d_fast
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only : sll_p_pi

  use sll_m_lagrange_interpolation_1d_fast, only : &
    sll_s_lagrange_interpolation_1d_fast_disp_fixed_no_bc, &
    sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodic, &
    sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodicl, &
    sll_s_lagrange_interpolation_1d_fast_disp_centered_periodicl


  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32 :: i, num_points
  sll_real64 :: alpha, xmin, xmax, l
  sll_real64, dimension(:), allocatable :: xi, fi, xp
  logical :: ctest

  ctest = .true.

  num_points = 100
  alpha = 0.2_f64
  ! known function values
  allocate(xi(1:num_points+1))
  allocate(fi(1:num_points+1))
  ! x values with offset
  allocate(xp(1:num_points+1))

  ! data initialization
  xmin = 0.0_f64
  xmax = num_points-1.0_f64
  l = xmax - xmin
  do i = 1, num_points+1
    xi(i) = real(i - 1, f64)
    fi(i) = f(xi(i), num_points)
    xp(i) = xi(i) + alpha
  end do

  
  call test_interpolation(num_points, fi(1:num_points), alpha, xp(1:num_points), 5, 0, 1.d-8, ctest)
  call test_interpolation(num_points, fi(1:num_points), alpha, xp(1:num_points), 3, 1, 8.d-6, ctest)
  call test_interpolation(num_points+1, fi, alpha, xp, 5, 2, 7.d-9, ctest)
  call test_interpolation(num_points+1, fi, alpha, xp, 4, 3, 3.d-7, ctest)
  call test_interpolation(num_points+1, fi, alpha, xp, 6, 3, 2.d-10, ctest)

  if ( ctest .eqv. .true. ) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if
  

  deallocate (xi,fi,xp)

contains

  subroutine test_interpolation(num_points, fi, alpha, xp, order, type, tolerance, ctest)
    sll_int32,  intent( in    ) :: num_points
    sll_real64, intent( in    ) :: fi(:)
    sll_real64, intent( in    ) :: alpha
    sll_real64, intent( in    ) :: xp(:)
    sll_int32,  intent( in    ) :: order
    sll_int32,  intent( in    ) :: type
    sll_real64, intent( in    ) :: tolerance
    logical,    intent( inout ) :: ctest

    character(len=2) :: pmessage
    sll_int32   :: num_cells
    sll_real64 :: fp(1:num_points)
    sll_real64 :: diff

    fp(:) = 0.0_f64
    diff = 0.0_f64

    if (type>1) then
       num_cells = num_points-1
    else
       num_cells = num_points
    end if

    write(pmessage, '(i2)') order

    if (type == 0) then ! no bc
       print*, "Test fixed_no_bc with order", pmessage, '.'
       call sll_s_lagrange_interpolation_1d_fast_disp_fixed_no_bc(fi, fp, alpha, order)
    elseif (type == 1) then ! periodic
       print*, "Test fixed_periodic with order", pmessage, '.'
       call sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodic(fi, fp, alpha, order)
    elseif (type == 2) then ! periodic with last value
       print*, "Test fixed_periodic_last with order", pmessage, '.'
       call sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodicl(fi, fp, alpha, order)
    elseif (type == 3) then ! periodic centered
       print*, "Test centered_periodic_last with order", pmessage, '.'
       call sll_s_lagrange_interpolation_1d_fast_disp_centered_periodicl(fi, fp, alpha, order)
    else
       print*, 'Interpolation type not implemented.'
    end if
       
    do i=1,num_points
       diff=max(diff, abs(f(xp(i),num_cells) - fp(i)))
    end do

    print*, "error =", diff

    if (diff > tolerance) then
       ctest = .false.
    end if

  end subroutine test_interpolation

  function f(x, num_points)
    sll_int32, intent(in) :: num_points
    sll_real64, intent(in) :: x
    sll_real64 :: f
    f = cos(2*sll_p_pi*x/num_points)
  end function

end program test_lagrange_interpolation_1d_fast
