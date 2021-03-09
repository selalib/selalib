program test_splines_pp_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_low_level_bsplines, only: &
      sll_s_uniform_bsplines_eval_basis

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_constants, only: &
      sll_p_twopi

   use sll_m_splines_pp, only: &
      sll_t_spline_pp_1d, &
      sll_s_spline_pp_init_1d, &
      sll_s_spline_pp_free_1d, &
      sll_s_spline_pp_b_to_pp_1d, &
      sll_f_spline_pp_horner_1d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_real64, parameter :: inv_2 = 1._f64/2._f64
   sll_real64, parameter :: inv_3 = 1._f64/3._f64
   sll_real64, parameter :: inv_4 = 1._f64/4._f64
   sll_real64, parameter :: inv_6 = 1._f64/6._f64
   sll_real64, parameter :: inv_8 = 1._f64/8._f64
   sll_real64, parameter :: inv_10 = 1._f64/10._f64
   sll_real64, parameter :: inv_24 = 1._f64/24._f64

   type(sll_t_spline_pp_1d) :: spline_pp
   sll_int32 :: degree
   sll_int32 :: n_cells
   logical   :: fail

   n_cells = 8
   do degree = 1, 5
      call spline_test(spline_pp, degree, n_cells, fail)

      if (fail .eqv. .false.) then
         write (*, *) 'PASSED'
      else
         write (*, *) 'FAILED'
         stop
      end if
   end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine spline_test(spline_pp, degree, n_cells, fail)
      type(sll_t_spline_pp_1d), intent(inout) :: spline_pp !arbitrary degree spline
      sll_int32, intent(in) :: degree !spline degree
      sll_int32, intent(in) :: n_cells !grid cells
      logical, intent(out)    :: fail

      sll_real64 :: b_coeffs(n_cells)
      sll_real64 :: pp_coeffs(degree + 1, n_cells)
      sll_real64 :: xp !particle position
      sll_real64 :: res
      sll_real64 :: res2
      sll_real64 :: domain(2) !left and right bound of the 1d grid
      sll_real64 :: delta_x !size of gridcells
      sll_real64, allocatable :: val(:)
      sll_int32  :: index1d
      sll_int32  :: index
      sll_real64 :: xi
      sll_int32  :: i

      integer, dimension(33) :: seed = [-1265139360, 1818735519, -1278717687, &
                                        -28571184, -2049390160, 2074167660, &
                                        -1778129250, -1663924455, -300142776, &
                                        1497205713, 1463052918, -1650171289, &
                                        1313784976, -1898838479, 2125570893, &
                                        -162457092, 1760636990, 524383974, &
                                        296008199, -171091367, 399322358, &
                                        967084750, 1776047718, -895222581, &
                                        -2070137937, -1280788435, 2086980348, &
                                        1463273178, 465948978, -701015021, &
                                        1313707185, 1192973109, 0]

      fail = .false.
      allocate (val(degree + 1))
      domain(1) = 0._f64
      domain(2) = sll_p_twopi
      delta_x = (domain(2) - domain(1))/real(n_cells, f64)

      call random_seed(put=seed)
      call random_number(b_coeffs)

      call sll_s_spline_pp_init_1d(spline_pp, degree, n_cells)

      call sll_s_spline_pp_b_to_pp_1d(spline_pp, n_cells, b_coeffs, pp_coeffs)

      call random_number(xp)
      xp = xp*(domain(2) - domain(1))

      ! calculate index of gridcell
      xi = (xp - domain(1))/delta_x
      index = floor(xi) + 1
      xi = xi - real(index - 1, f64)
      res = sll_f_spline_pp_horner_1d(degree, pp_coeffs(:, :), xi, index)

      index = index - degree

      call sll_s_uniform_bsplines_eval_basis(degree, xi, val)

      res2 = 0.0_f64
      do i = 1, degree + 1
         index1d = modulo(index + i - 2, n_cells) + 1
         res2 = res2 + b_coeffs(index1d)*val(i)
      end do
      !write(*,*) 'Fehler horner vs normal:', abs(res-res2)
      if (abs(res - res2) > 1d-15) fail = .true.
      if (real(spline_pp%degree - degree) > 1E-15) then
         fail = .true.
         print *, 'error in evaluate'
      end if

      !test horner for arbitrary polynomials
      call random_number(xp)
      res = sll_f_spline_pp_horner_1d(degree, pp_coeffs, xp, 1)
      res2 = 0._f64
      do i = 1, degree + 1
         res2 = res2 + pp_coeffs(i, 1)*xp**((degree + 1) - i)
      end do
      if (abs(res - res2) > 1d-12) then
         fail = .true.
         print *, xp
         print *, 'error in horner'
      end if
      ! print*, 'res-res2=', res-res2

      deallocate (val)
      call sll_s_spline_pp_free_1d(spline_pp)

   end subroutine spline_test

end program test_splines_pp_1d
