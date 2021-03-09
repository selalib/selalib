program test_dg_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_dg_interpolator_1d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32, parameter  :: n_cells = 32
   sll_int32, parameter  :: degree = 3
   sll_int32, parameter  :: n = n_cells*degree
   sll_int32, parameter  :: m = 512
   sll_real64, parameter :: tol = 6.d-5 ! tolerance for error
   sll_real64            :: error1

   call test_disp_for_given_n(error1)

contains

   subroutine test_disp_for_given_n(error)
      sll_real64, intent(out) :: error

      type(sll_t_dg_interpolator_1d) :: interp

      sll_real64, allocatable, dimension(:) :: point
      sll_real64, allocatable, dimension(:) :: pdata
      sll_real64, allocatable, dimension(:) :: fdata
      sll_real64, allocatable, dimension(:) :: phdata
      sll_real64, allocatable, dimension(:) :: gdata
      sll_real64, allocatable, dimension(:, :) :: weights1
      sll_real64, allocatable, dimension(:, :) :: weights2

      sll_int32 :: ierr, i

      sll_real64  :: x_min, x_max, delta, alpha
      SLL_ALLOCATE(pdata(n), ierr)
      SLL_ALLOCATE(phdata(n + degree), ierr)
      SLL_ALLOCATE(fdata(n), ierr)
      SLL_ALLOCATE(gdata(n), ierr)
      SLL_ALLOCATE(weights1(degree, degree), ierr)
      SLL_ALLOCATE(weights2(degree, degree), ierr)

      print *, 'Initialize data and point array'
      x_min = 0.0_f64
      x_max = 2.0_f64*sll_p_pi

      print *, 'DG interpolation'
      call sll_s_dg_interpolator_1d_init &
         (interp, n_cells, degree, x_min, x_max, sll_p_dg_gauss_legendre)

      delta = (x_max - x_min)/real(n_cells, f64)

      alpha = -delta*1.2_f64

      do i = 1, n
         pdata(i) = f(interp%positions(i))
         gdata(i) = f(interp%positions(i) + alpha)
      end do

      phdata(degree*2 + 1:degree + n) = pdata(1:n - degree)
      phdata(1:degree*2) = pdata(n - 2*degree + 1:n)

      call sll_s_dg_interpolator_1d_interpolate_array_disp_periodic(interp, n, pdata, alpha, fdata)

      error = maxval(abs(gdata - fdata))
      print *, 'error=', error

      call sll_s_dg_interpolator_1d_interpolate_array_disp_halo(interp, phdata, 0.8_f64, weights1, weights2, fdata)

      error = maxval(abs(gdata - fdata))
      print *, 'error=', error

      if (error < tol) then
         print *, "PASSED."
      else
         print *, "FAILED."
      end if

   end subroutine test_disp_for_given_n

   function f(x)

      sll_real64 :: x
      sll_real64 :: f

      f = 2.0_f64*(sin(x) + 2.5_f64 + cos(x))

   end function f

end program test_dg_interpolator_1d
