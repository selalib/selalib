!> Test programm for the splines_pp in 2d (general function)
!> author: Katharina Kormann, IPP

program test_splines_pp_2d_boundary
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   ! use sll_m_low_level_bsplines, only: &
   !      sll_s_uniform_bsplines_eval_basis

   use sll_m_bsplines_non_uniform, only: &
      sll_t_bsplines_non_uniform

   use sll_m_constants, only: &
      sll_p_twopi

   use sll_m_splines_pp, only: &
      sll_t_spline_pp_2d, &
      sll_s_spline_pp_init_2d, &
      sll_s_spline_pp_free_2d, &
      sll_s_spline_pp_b_to_pp_2d, &
      sll_f_spline_pp_horner_2d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_real64, parameter :: inv_2 = 1._f64/2._f64
   sll_real64, parameter :: inv_3 = 1._f64/3._f64
   sll_real64, parameter :: inv_4 = 1._f64/4._f64
   sll_real64, parameter :: inv_6 = 1._f64/6._f64
   sll_real64, parameter :: inv_8 = 1._f64/8._f64
   sll_real64, parameter :: inv_10 = 1._f64/10._f64
   sll_real64, parameter :: inv_24 = 1._f64/24._f64

   type(sll_t_spline_pp_2d) :: spline_pp
   sll_int32 :: degree(2)
   sll_int32 :: n_cells(2)
   sll_int32 :: boundary(2)
   sll_real64 :: xnorm(2)
   logical   :: fail

   sll_int32 :: nseed
   sll_int32, allocatable :: seed(:)

   call random_seed(size=nseed)
   allocate (seed(nseed))
   seed = 42
   call random_seed(put=seed)

   fail = .false.

   n_cells = [50, 50]
   call random_number(xnorm)
   print *, 'Order 3'
   degree = [3, 3]
   print *, 'periodic'
   boundary = [0, 0]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'clamped, first interval, second last interval'
   boundary = [1, 1]
   xnorm = [0.01_f64, 0.97_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'clamped, last interval, second interval'
   xnorm = [0.99_f64, 0.03_f64]
   boundary = [1, 1]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'clamped, periodic'
   boundary = [1, 0]
   xnorm = [0.035_f64, 0.567_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'clamped, periodic'
   boundary = [1, 0]
   xnorm = [0.567_f64, 0.035_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'clamped-clampeddiri, periodic'
   degree = [3, 3]
   boundary = [2, 0]
   xnorm = [0.985_f64, 0.567_f64]
   print *, 'clampeddiri, periodic, first interval'
   degree = [3, 3]
   boundary = [3, 0]
   xnorm = [0.085_f64, 0.567_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'clampeddiri, periodic, last interval'
   degree = [3, 3]
   boundary = [3, 0]
   xnorm = [0.985_f64, 0.567_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'clampeddiri-clamped, periodic'
   degree = [3, 3]
   boundary = [4, 0]
   xnorm = [0.085_f64, 0.567_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'periodic, clamped-clampeddiri'
   degree = [3, 3]
   boundary = [0, 2]
   xnorm = [0.567_f64, 0.995_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'periodic, clampeddiri, last interval'
   degree = [3, 3]
   boundary = [0, 3]
   xnorm = [0.567_f64, 0.99_f64]!0.975_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'periodic, clampeddiri, first interval'
   degree = [3, 3]
   boundary = [0, 3]
   xnorm = [0.567_f64, 0.01_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'periodic, clampeddiri-clamped, last interval'
   degree = [3, 3]
   boundary = [0, 4]
   xnorm = [0.567_f64, 0.995_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'
   print *, 'periodic, clampeddiri-clamped, first interval'
   degree = [3, 3]
   boundary = [0, 4]
   xnorm = [0.567_f64, 0.01_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
   if (fail) stop 'FAILED'

   print *, 'Order 2'
   degree = [2, 2]
   print *, 'periodic, clamped'
   boundary = [0, 1]
   xnorm = [0.767_f64, 0.99_f64]
   call spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)

   if (fail) stop 'FAILED'

   write (*, *) 'PASSED'

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine spline_test(spline_pp, degree, n_cells, boundary, xnorm, fail)
      type(sll_t_spline_pp_2d), intent(inout) :: spline_pp !arbitrary degree spline
      sll_int32, intent(in) :: degree(2) !spline degree
      sll_int32, intent(in) :: n_cells(2) !grid cells
      sll_int32, intent(in) :: boundary(2) ! boundaries
      sll_real64, intent(in) :: xnorm(2)
      logical, intent(out)    :: fail

      sll_real64 :: b_coeffs((n_cells(1) + degree(1))*(n_cells(2) + degree(2)))
      sll_real64 :: b_coeffs2d((n_cells(1) + degree(1)), (n_cells(2) + degree(2)))
      sll_real64 :: b_coeffs_cp((n_cells(1) + degree(1))*(n_cells(2) + degree(2)))
      sll_real64 :: pp_coeffs((degree(1) + 1)*(degree(2) + 1), n_cells(1)*n_cells(2))
      sll_real64 :: xp(2) !particle position
      sll_real64 :: res
      sll_real64 :: res2
      sll_real64 :: domain(2, 2) !left and right bound of the 1d grid
      sll_real64 :: delta_x(2) !size of gridcells
      sll_real64, allocatable :: val1(:)
      sll_real64, allocatable :: val2(:)
      sll_int32 :: index2d
      sll_int32 :: index1d(2)
      sll_int32 :: indices(2)
      sll_real64 :: xi(2)
      sll_int32 :: i, j, ind1, ind2, n_coeffs(2), ind, low, lowm
      type(sll_t_bsplines_non_uniform) :: bspline1
      type(sll_t_bsplines_non_uniform) :: bspline2
      sll_real64, allocatable :: breaks1(:)
      sll_real64, allocatable :: breaks2(:)

      !fail=.false.
      allocate (val1(degree(1) + 1))
      allocate (val2(degree(2) + 1))
      domain(1, :) = 0._f64
      domain(2, :) = sll_p_twopi
      delta_x = (domain(2, :) - domain(1, :))/real(n_cells, f64)

      allocate (breaks1(n_cells(1) + 1))
      breaks1(1) = domain(1, 1)
      do i = 1, n_cells(1)
         breaks1(i + 1) = breaks1(i) + delta_x(1)
      end do
      allocate (breaks2(n_cells(2) + 1))
      breaks2(1) = domain(1, 2)
      do i = 1, n_cells(2)
         breaks2(i + 1) = breaks2(i) + delta_x(2)
      end do

      call random_seed(put=seed)
      call random_number(b_coeffs)

      call sll_s_spline_pp_init_2d(spline_pp, degree, n_cells, boundary)
      call sll_s_spline_pp_b_to_pp_2d(spline_pp, n_cells, b_coeffs, pp_coeffs)

      n_coeffs(1) = spline_pp%spline1%n_coeffs
      n_coeffs(2) = spline_pp%spline2%n_coeffs
      !print*, n_coeffs

      !write(34,*) b_coeffs
      !print*, b_coeffs(1), b_coeffs(1+n_coeffs(1))

      !write(32,*) pp_coeffs
      !print*, pp_coeffs(:,1)

      b_coeffs_cp = b_coeffs

      ind = 1
      do j = 1, n_cells(2) + degree(2)
         do i = 1, n_cells(1) + degree(1)
            b_coeffs2d(i, j) = b_coeffs(ind)
            ind = ind + 1
         end do
      end do

      if (boundary(1) == 2 .and. boundary(2) < 2) then
         low = 0
         lowm = 0
         do i = 1, n_coeffs(2)
            b_coeffs_cp(low + 1:low + n_coeffs(1)) = b_coeffs(lowm + 1:lowm + n_coeffs(1))
            b_coeffs_cp(low + n_coeffs(1) + 1) = 0.0_f64
            low = low + n_coeffs(1) + 1
            lowm = lowm + n_coeffs(1)
         end do

         n_coeffs(1) = n_coeffs(1) + 1
      elseif (boundary(1) == 3 .and. boundary(2) < 2) then
         low = 1
         lowm = 0
         do i = 1, n_coeffs(2)
            b_coeffs_cp(low) = 0.0_f64
            b_coeffs_cp(low + 1:low + n_coeffs(1)) = b_coeffs(lowm + 1:lowm + n_coeffs(1))
            b_coeffs_cp(low + n_coeffs(1) + 1) = 0.0_f64
            low = low + n_coeffs(1) + 2
            lowm = lowm + n_coeffs(1)
         end do

         n_coeffs(1) = n_coeffs(1) + 2
      elseif (boundary(1) == 4 .and. boundary(2) < 2) then
         low = 1
         lowm = 0
         do i = 1, n_coeffs(2)
            b_coeffs_cp(low) = 0.0_f64
            b_coeffs_cp(low + 1:low + n_coeffs(1)) = b_coeffs(lowm + 1:lowm + n_coeffs(1))
            low = low + n_coeffs(1) + 1
            lowm = lowm + n_coeffs(1)
         end do

         n_coeffs(1) = n_coeffs(1) + 1

      elseif (boundary(2) == 2 .and. boundary(1) < 2) then
         b_coeffs_cp(n_coeffs(2)*n_coeffs(1) + 1:(n_coeffs(2) + 1)*n_coeffs(1)) = 0.0_f64

         n_coeffs(2) = n_coeffs(2) + 1
      elseif (boundary(2) == 3 .and. boundary(1) < 2) then
         b_coeffs_cp(1:n_coeffs(1)) = 0.0_f64
         b_coeffs_cp(n_coeffs(1) + 1:(n_coeffs(2) + 1)*n_coeffs(1)) = b_coeffs(1:n_coeffs(2)*n_coeffs(1))
         b_coeffs_cp((n_coeffs(2) + 1)*n_coeffs(1) + 1:(n_coeffs(2) + 2)*n_coeffs(1)) = 0.0_f64

         n_coeffs(2) = n_coeffs(2) + 2
      elseif (boundary(2) == 4 .and. boundary(1) < 2) then
         b_coeffs_cp(1:n_coeffs(1)) = 0.0_f64
         b_coeffs_cp(n_coeffs(1) + 1:(n_coeffs(2) + 1)*n_coeffs(1)) = b_coeffs(1:n_coeffs(2)*n_coeffs(1))

         n_coeffs(2) = n_coeffs(2) + 1
      end if

      !print*, 'ncoef', n_coeffs

      ! Initialize the bspline
      if (boundary(1) == 0) then
         call bspline1%init(degree(1), .true., breaks1)
      else
         call bspline1%init(degree(1), .false., breaks1)
      end if
      if (boundary(2) == 0) then
         call bspline2%init(degree(2), .true., breaks2)
      else
         call bspline2%init(degree(2), .false., breaks2)
      end if

      !call random_seed()
      !call random_number(xp)
      !xp=[0.8_f64,0.2_f64]
      xp = xnorm*(domain(2, :) - domain(1, :))
      !xp(1) = 0.1_f64
      !xp(2) = 6.2_f64
      !xp(1) = xi(2)
      !xp(2) = xi(1)

      ! calculate index of gridcell
      xi = (xp - domain(1, :))/delta_x
      indices = floor(xi) + 1
      xi = xi - real(indices - 1, f64)
      res = sll_f_spline_pp_horner_2d(degree, pp_coeffs, xi, indices, n_cells)

      if (boundary(1) == 0) then
         indices(1) = indices(1) - degree(1)
      end if
      if (boundary(2) == 0) then
         indices(2) = indices(2) - degree(2)
      end if

      call bspline1%eval_basis(xp(1), val1(:), ind1)
      !print*, breaks1
      !print*, xp

      call bspline2%eval_basis(xp(2), val2(:), ind2)
      !call sll_s_uniform_bsplines_eval_basis(degree(1), xi(1), val1(:))
      !call sll_s_uniform_bsplines_eval_basis(degree(2), xi(2), val2(:))

      ! print*, indices, 'ind'
      res2 = 0.0_f64
      do i = 1, degree(1) + 1
         if (boundary(1) == 0) then
            index1d(1) = modulo(indices(1) + i - 2, n_cells(1)) + 1
         else
            index1d(1) = indices(1) + i - 1
         end if

         do j = 1, degree(2) + 1
            if (boundary(2) == 0) then
               index1d(2) = modulo(indices(2) + j - 2, n_cells(2))
            else
               index1d(2) = indices(2) + j - 2
            end if
            index2d = index1d(1) + index1d(2)*n_coeffs(1)
            !print*, index1d, n_coeffs(1)
            !print*, index2d, b_coeffs_cp(index2d)
            res2 = res2 + b_coeffs_cp(index2d)*val1(i)*val2(j)
         end do
      end do
      !print*, res, res2
      !print*, 'res-res2=', res-res2
      !write(*,*) 'Fehler horner vs normal:', abs(res-res2)
      if (abs(res - res2) > 5E-14) then
         fail = .true.
         print *, 'error in evaluate', res - res2
      end if

      !test horner for arbitrary polynomials
      call random_seed(put=seed)
      call random_number(xp)
      res = sll_f_spline_pp_horner_2d(degree, pp_coeffs, xp, [1, 1], [1, 1])
      res2 = 0._f64
      do i = 1, degree(1) + 1
         do j = 1, degree(2) + 1
            res2 = res2 + pp_coeffs(i + (j - 1)*(degree(1) + 1), 1)*xp(1)**((degree(1) + 1) - i)*xp(2)**((degree(2) + 1) - j)
         end do
      end do
      if (abs(res - res2) > 1E-12) then
         fail = .true.
         print *, xp
         print *, 'error in horner'
      end if

      !print*, 'res-res2=', res-res2

      deallocate (val1)
      deallocate (val2)
      call sll_s_spline_pp_free_2d(spline_pp)

   end subroutine spline_test

end program test_splines_pp_2d_boundary
