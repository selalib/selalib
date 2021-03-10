!> @brief
!> module for evaluation of a polynomial
!> in pp form using Horner's algorithm
!
!> @authors  Cèline Caldini-Queiros
!> @authors  Yaman Güçlü, IPP Garching
!
module sll_m_horner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

   use sll_m_working_precision, only: f64

   implicit none

   public :: &
      sll_f_horner_1d_eval, &
      sll_f_horner_2d_eval, &
      sll_f_horner_3d_eval

   private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   !> @brief       Horner's algorithm in 1D
   !> @description evaluation of a polynomial or its derivative
   !>
   !> @param[in]  coeffs: pp coeffs
   !> @param[in]  x     : point of evaluation
   !> @param[in]  deriv : order of derivative
   !> @param[out] res   : value
   !-----------------------------------------------------------------------------
   SLL_PURE function sll_f_horner_1d_eval(coeffs, x, deriv) result(res)
      real(wp), intent(in) :: coeffs(:)
      real(wp), intent(in) :: x
      integer, intent(in) :: deriv
      real(wp) :: res

      real(wp) :: fmmjdr
      integer  :: i1, n1

      n1 = ubound(coeffs, 1)

      fmmjdr = real(n1 - 1 - deriv, wp)
      res = coeffs(n1)
      do i1 = 1, n1 - 1 - deriv
         res = (res/fmmjdr)*x + coeffs(n1 - i1)
         fmmjdr = fmmjdr - 1.0_wp
      end do

   end function sll_f_horner_1d_eval

   !-----------------------------------------------------------------------------
   !> @brief       Horner's algorithm in 2D
   !> @description evaluation of a polynomial or its (mixed) partial derivative
   !>
   !> @param[in]  coeffs: pp coeffs
   !> @param[in]  x     : point of evaluation          (2 components)
   !> @param[in]  deriv : order of partial derivatives (2 components)
   !> @param[out] res   : value
   !-----------------------------------------------------------------------------
   SLL_PURE function sll_f_horner_2d_eval(coeffs, x, deriv) result(res)
      real(wp), intent(in) :: coeffs(:, :)
      real(wp), intent(in) ::      x(:)
      integer, intent(in) ::  deriv(:)
      real(wp) :: res

      real(wp) :: coeffi(ubound(coeffs, 2))
      integer  :: i2, n2

      n2 = ubound(coeffs, 2)

      do i2 = 1, n2
         coeffi(i2) = sll_f_horner_1d_eval(coeffs(:, i2), x(1), deriv(1))
      end do

      res = sll_f_horner_1d_eval(coeffi, x(2), deriv(2))

   end function sll_f_horner_2d_eval

   !-----------------------------------------------------------------------------
   !> @brief       Horner's algorithm in 3D
   !> @description evaluation of a polynomial or its (mixed) partial derivative
   !>
   !> @param[in]  coeffs: pp coeffs
   !> @param[in]  x     : point of evaluation          (3 components)
   !> @param[in]  deriv : order of partial derivatives (3 components)
   !> @param[out] res   : value
   !-----------------------------------------------------------------------------
   SLL_PURE function sll_f_horner_3d_eval(coeffs, x, deriv) result(res)
      real(wp), intent(in) :: coeffs(:, :, :)
      real(wp), intent(in) ::      x(:)
      integer, intent(in) ::  deriv(:)
      real(wp) :: res

      real(wp) :: coeffi(ubound(coeffs, 3))
      integer  :: i3, n3

      n3 = ubound(coeffs, 3)

      do i3 = 1, n3
         coeffi(i3) = sll_f_horner_2d_eval(coeffs(:, :, i3), x(1:2), deriv(1:2))
      end do

      res = sll_f_horner_1d_eval(coeffi, x(3), deriv(3))

   end function sll_f_horner_3d_eval

end module sll_m_horner
