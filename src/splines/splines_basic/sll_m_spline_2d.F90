!> @ingroup splines
!> @brief   Module for tensor-product 2D splines.
!>
!> @details Here we define a 2D tensor-product spline type as an element of the
!>          linear space given by the span of 2D basis functions, which in turn
!>          are obtained as tensor-product of 1D B-splines.
!>          A 2D tensor-product B-splines type is not implemented, because an
!>          object of this type is completely described by the combination of
!>          two separate 1D B-splines objects, but this decision may be changed.
!>          Therefore, initialization of a 2D spline object requires two existing
!>          B-splines objects, to which two private (polymorphic) pointers are
!>          associated.
!>          The B-spline coefficients are stored in a 2D public allocatable array;
!>          at initialization the array is allocated to the proper shape and all
!>          values are set to zero.
!>          In most situations the B-spline coefficients are not set directly by
!>          the end user, but are computed by some other object (e.g., a Poisson
!>          solver or a spline interpolator).
!>          Various public methods allow the user to evaluate the 2D spline
!>          S(x1,x2) and its partial derivatives ∂S(x1,x2)/∂x1 and ∂S(x1,x2)/∂x2
!>          at any position (x1,x2).
!>
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

   use sll_m_working_precision, only: f64
   use sll_m_bsplines_base, only: sll_c_bsplines

   implicit none

   public :: &
      sll_t_spline_2d

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !> 2D tensor-product spline
   type :: sll_t_spline_2d

      real(wp), allocatable      :: bcoef(:, :)
      class(sll_c_bsplines), pointer, private :: bspl1 => null()
      class(sll_c_bsplines), pointer, private :: bspl2 => null()

   contains

      procedure :: init => s_spline_2d__init
      procedure :: free => s_spline_2d__free
      procedure :: belongs_to_space => f_spline_2d__belongs_to_space
      procedure :: eval => f_spline_2d__eval
      procedure :: eval_deriv_x1 => f_spline_2d__eval_deriv_x1
      procedure :: eval_deriv_x2 => f_spline_2d__eval_deriv_x2
      procedure :: eval_deriv_x1x2 => f_spline_2d__eval_deriv_x1x2
      procedure :: eval_array => s_spline_2d__eval_array
      procedure :: eval_array_deriv_x1 => s_spline_2d__eval_array_deriv_x1
      procedure :: eval_array_deriv_x2 => s_spline_2d__eval_array_deriv_x2
      procedure :: eval_array_deriv_x1x2 => s_spline_2d__eval_array_deriv_x1x2

   end type sll_t_spline_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   !> @brief      Initialize 2D spline object as element of span(B-splines)
   !> @param[out] self         2D spline: new element of tensor-product space
   !> @param[in]  bsplines_x1  B-splines along x1
   !> @param[in]  bsplines_x2  B-splines along x2
   !-----------------------------------------------------------------------------
   subroutine s_spline_2d__init(self, bsplines_x1, bsplines_x2)

      class(sll_t_spline_2d), intent(out)         :: self
      class(sll_c_bsplines), intent(in), target :: bsplines_x1
      class(sll_c_bsplines), intent(in), target :: bsplines_x2

      ! Store pointer to B-splines
      self%bspl1 => bsplines_x1
      self%bspl2 => bsplines_x2

      ! Allocate array of spline coefficients: in case of periodic BCs, the last
      ! p coefficients are a periodic copy of the first p ones (p=p1,p2)
      associate (n1 => bsplines_x1%ncells, &
                 n2 => bsplines_x2%ncells, &
                 p1 => bsplines_x1%degree, &
                 p2 => bsplines_x2%degree)

         allocate (self%bcoef(1:n1 + p1, 1:n2 + p2))

      end associate

      ! Set all coefficients to zero
      self%bcoef = 0.0_f64

   end subroutine s_spline_2d__init

   !-----------------------------------------------------------------------------
   !> @brief        Destroy 2D spline (re-initialization is possible afterwards)
   !> @param[inout] self  2D spline
   !-----------------------------------------------------------------------------
   subroutine s_spline_2d__free(self)

      class(sll_t_spline_2d), intent(inout) :: self

      deallocate (self%bcoef)
      nullify (self%bspl1)
      nullify (self%bspl2)

   end subroutine s_spline_2d__free

   !-----------------------------------------------------------------------------
   !> @brief     Check if 2D spline belongs to tensor-product span of B-splines
   !> @param[in] self         2D spline
   !> @param[in] bsplines_x1  B-splines along x1
   !> @param[in] bsplines_x2  B-splines along x2
   !> @returns                true/false answer
   !-----------------------------------------------------------------------------
   SLL_PURE function f_spline_2d__belongs_to_space(self, bsplines_x1, bsplines_x2) result(in_space)

      class(sll_t_spline_2d), intent(in)         :: self
      class(sll_c_bsplines), intent(in), target :: bsplines_x1
      class(sll_c_bsplines), intent(in), target :: bsplines_x2
      logical :: in_space

      in_space = associated(self%bspl1, bsplines_x1) .and. &
                 associated(self%bspl2, bsplines_x2)

   end function f_spline_2d__belongs_to_space

   !-----------------------------------------------------------------------------
   !> @brief     Evaluate value of 2D spline at location (x1,x2)
   !> @param[in] self  2D spline
   !> @param[in] x1    x1 coordinate of evaluation point
   !> @param[in] x2    x2 coordinate of evaluation point
   !> @returns         spline value y=S(x1,x2)
   !-----------------------------------------------------------------------------
   SLL_PURE function f_spline_2d__eval(self, x1, x2) result(y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1
      real(wp), intent(in) :: x2
      real(wp) :: y

      integer :: jmin(2)
      integer :: jmax(2)
      integer :: k1, k2

      ! Automatic arrays
      real(wp) :: values1(1:self%bspl1%degree + 1)
      real(wp) :: values2(1:self%bspl2%degree + 1)

      ! Compute arrays v1 and v2 of B-spline values
      call self%bspl1%eval_basis(x1, values1, jmin(1))
      call self%bspl2%eval_basis(x2, values2, jmin(2))

      jmax(1) = jmin(1) + self%bspl1%degree
      jmax(2) = jmin(2) + self%bspl2%degree

      ! Determine (matrix) block C of B-spline coefficients to be used
      ! and compute scalar product <v1,v2> = (v1^T)*C*v2
      associate (bcoef => self%bcoef(jmin(1):jmax(1), jmin(2):jmax(2)))
         y = 0.0_f64
         do k2 = 1, 1 + self%bspl2%degree
            do k1 = 1, 1 + self%bspl1%degree
               y = y + bcoef(k1, k2)*values1(k1)*values2(k2)
            end do
         end do
      end associate

   end function f_spline_2d__eval

   !-----------------------------------------------------------------------------
   !> @brief     Evaluate x1-derivative of 2D spline at location (x1,x2)
   !> @param[in] self  2D spline
   !> @param[in] x1    x1 coordinate of evaluation point
   !> @param[in] x2    x2 coordinate of evaluation point
   !> @returns         spline partial derivative y=∂S(x1,x2)/∂x1
   !-----------------------------------------------------------------------------
   SLL_PURE function f_spline_2d__eval_deriv_x1(self, x1, x2) result(y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1
      real(wp), intent(in) :: x2
      real(wp) :: y

      integer :: jmin(2)
      integer :: jmax(2)
      integer :: k1, k2

      ! Automatic arrays
      real(wp) :: derivs1(1:self%bspl1%degree + 1)
      real(wp) :: values2(1:self%bspl2%degree + 1)

      ! Compute arrays d1 and v2 of B-spline derivatives/values
      call self%bspl1%eval_deriv(x1, derivs1, jmin(1))
      call self%bspl2%eval_basis(x2, values2, jmin(2))

      jmax(1) = jmin(1) + self%bspl1%degree
      jmax(2) = jmin(2) + self%bspl2%degree

      ! Determine (matrix) block C of B-spline coefficients to be used
      ! and compute scalar product <d1,v2> = (d1^T)*C*v2
      associate (bcoef => self%bcoef(jmin(1):jmax(1), jmin(2):jmax(2)))
         y = 0.0_f64
         do k2 = 1, 1 + self%bspl2%degree
            do k1 = 1, 1 + self%bspl1%degree
               y = y + bcoef(k1, k2)*derivs1(k1)*values2(k2)
            end do
         end do
      end associate

   end function f_spline_2d__eval_deriv_x1

   !-----------------------------------------------------------------------------
   !> @brief     Evaluate x2-derivative of 2D spline at location (x1,x2)
   !> @param[in] self  2D spline
   !> @param[in] x1    x1 coordinate of evaluation point
   !> @param[in] x2    x2 coordinate of evaluation point
   !> @returns         spline partial derivative y=∂S(x1,x2)/∂x2
   !-----------------------------------------------------------------------------
   SLL_PURE function f_spline_2d__eval_deriv_x2(self, x1, x2) result(y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1
      real(wp), intent(in) :: x2
      real(wp) :: y

      integer :: jmin(2)
      integer :: jmax(2)
      integer :: k1, k2

      ! Automatic arrays
      real(wp) :: values1(1:self%bspl1%degree + 1)
      real(wp) :: derivs2(1:self%bspl2%degree + 1)

      ! Compute arrays v1 and d2 of B-spline values/derivatives
      call self%bspl1%eval_basis(x1, values1, jmin(1))
      call self%bspl2%eval_deriv(x2, derivs2, jmin(2))

      jmax(1) = jmin(1) + self%bspl1%degree
      jmax(2) = jmin(2) + self%bspl2%degree

      ! Determine (matrix) block C of B-spline coefficients to be used
      ! and compute scalar product <v1,d2> = (v1^T)*C*d2
      associate (bcoef => self%bcoef(jmin(1):jmax(1), jmin(2):jmax(2)))
         y = 0.0_f64
         do k2 = 1, 1 + self%bspl2%degree
            do k1 = 1, 1 + self%bspl1%degree
               y = y + bcoef(k1, k2)*values1(k1)*derivs2(k2)
            end do
         end do
      end associate

   end function f_spline_2d__eval_deriv_x2

   !-----------------------------------------------------------------------------
   !> @brief     Evaluate x1-x2 mixed derivative of 2D spline at location (x1,x2)
   !> @param[in] self  2D spline
   !> @param[in] x1    x1 coordinate of evaluation point
   !> @param[in] x2    x2 coordinate of evaluation point
   !> @returns         spline partial derivative y = ∂^2/(∂x1 ∂x2) S(x1,x2)
   !-----------------------------------------------------------------------------
   SLL_PURE function f_spline_2d__eval_deriv_x1x2(self, x1, x2) result(y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1
      real(wp), intent(in) :: x2
      real(wp) :: y

      integer :: jmin(2)
      integer :: jmax(2)
      integer :: k1, k2

      ! Automatic arrays
      real(wp) :: derivs1(1:self%bspl1%degree + 1)
      real(wp) :: derivs2(1:self%bspl2%degree + 1)

      ! Compute arrays d1 and d2 of B-spline values/derivatives
      call self%bspl1%eval_deriv(x1, derivs1, jmin(1))
      call self%bspl2%eval_deriv(x2, derivs2, jmin(2))

      jmax(1) = jmin(1) + self%bspl1%degree
      jmax(2) = jmin(2) + self%bspl2%degree

      ! Determine (matrix) block C of B-spline coefficients to be used
      ! and compute scalar product <d1,d2> = (d1^T)*C*d2
      associate (bcoef => self%bcoef(jmin(1):jmax(1), jmin(2):jmax(2)))
         y = 0.0_f64
         do k2 = 1, 1 + self%bspl2%degree
            do k1 = 1, 1 + self%bspl1%degree
               y = y + bcoef(k1, k2)*derivs1(k1)*derivs2(k2)
            end do
         end do
      end associate

   end function f_spline_2d__eval_deriv_x1x2

   !-----------------------------------------------------------------------------
   !> @brief      Evaluate value of 2D spline at multiple locations
   !> @param[in]  self  2D spline
   !> @param[in]  x1    2D array of x1 coordinates of evaluation points
   !> @param[in]  x2    2D array of x2 coordinates of evaluation points
   !> @param[out] y     2D array of spline values y[i,j]=S(x1[i,j],x2[i,j])
   !-----------------------------------------------------------------------------
   SLL_PURE subroutine s_spline_2d__eval_array(self, x1, x2, y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1(:, :)
      real(wp), intent(in) :: x2(:, :)
      real(wp), intent(out) :: y(:, :)

      integer :: i1, i2

      SLL_ASSERT(size(x1, 1) == size(y, 1))
      SLL_ASSERT(size(x1, 2) == size(y, 2))
      SLL_ASSERT(size(x2, 1) == size(y, 1))
      SLL_ASSERT(size(x2, 2) == size(y, 2))

      do i2 = 1, size(x2, 2)
         do i1 = 1, size(x1, 1)
            y(i1, i2) = f_spline_2d__eval(self, x1(i1, i2), x2(i1, i2))
         end do
      end do

   end subroutine s_spline_2d__eval_array

   !-----------------------------------------------------------------------------
   !> @brief      Evaluate x1-derivative of 2D spline at multiple locations
   !> @param[in]  self  2D spline
   !> @param[in]  x1    2D array of x1 coordinates of evaluation points
   !> @param[in]  x2    2D array of x2 coordinates of evaluation points
   !> @param[out] y     2D array of spline partial derivatives
   !>                   y[i,j]=∂S(x1[i,j],x2[i,j])/∂x1
   !-----------------------------------------------------------------------------
   SLL_PURE subroutine s_spline_2d__eval_array_deriv_x1(self, x1, x2, y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1(:, :)
      real(wp), intent(in) :: x2(:, :)
      real(wp), intent(out) :: y(:, :)

      integer :: i1, i2

      SLL_ASSERT(size(x1, 1) == size(y, 1))
      SLL_ASSERT(size(x1, 2) == size(y, 2))
      SLL_ASSERT(size(x2, 1) == size(y, 1))
      SLL_ASSERT(size(x2, 2) == size(y, 2))

      do i2 = 1, size(x2, 2)
         do i1 = 1, size(x1, 1)
            y(i1, i2) = f_spline_2d__eval_deriv_x1(self, x1(i1, i2), x2(i1, i2))
         end do
      end do

   end subroutine s_spline_2d__eval_array_deriv_x1

   !-----------------------------------------------------------------------------
   !> @brief      Evaluate x2-derivative of 2D spline at multiple locations
   !> @param[in]  self  2D spline
   !> @param[in]  x1    2D array of x1 coordinates of evaluation points
   !> @param[in]  x2    2D array of x2 coordinates of evaluation points
   !> @param[out] y     2D array of spline partial derivatives
   !>                   y[i,j]=∂S(x1[i,j],x2[i,j])/∂x2
   !-----------------------------------------------------------------------------
   SLL_PURE subroutine s_spline_2d__eval_array_deriv_x2(self, x1, x2, y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1(:, :)
      real(wp), intent(in) :: x2(:, :)
      real(wp), intent(out) :: y(:, :)

      integer :: i1, i2

      SLL_ASSERT(size(x1, 1) == size(y, 1))
      SLL_ASSERT(size(x1, 2) == size(y, 2))
      SLL_ASSERT(size(x2, 1) == size(y, 1))
      SLL_ASSERT(size(x2, 2) == size(y, 2))

      do i2 = 1, size(x2, 2)
         do i1 = 1, size(x1, 1)
            y(i1, i2) = f_spline_2d__eval_deriv_x2(self, x1(i1, i2), x2(i1, i2))
         end do
      end do

   end subroutine s_spline_2d__eval_array_deriv_x2

   !-----------------------------------------------------------------------------
   !> @brief      Evaluate x1-x2 mixed derivative of 2D spline at multiple locs.
   !> @param[in]  self  2D spline
   !> @param[in]  x1    2D array of x1 coordinates of evaluation points
   !> @param[in]  x2    2D array of x2 coordinates of evaluation points
   !> @param[out] y     2D array of spline partial derivatives
   !>                   y[i,j] = ∂^2/(∂x1 ∂x2) S(x1[i,j],x2[i,j])
   !-----------------------------------------------------------------------------
   SLL_PURE subroutine s_spline_2d__eval_array_deriv_x1x2(self, x1, x2, y)

      class(sll_t_spline_2d), intent(in) :: self
      real(wp), intent(in) :: x1(:, :)
      real(wp), intent(in) :: x2(:, :)
      real(wp), intent(out) :: y(:, :)

      integer :: i1, i2

      SLL_ASSERT(size(x1, 1) == size(y, 1))
      SLL_ASSERT(size(x1, 2) == size(y, 2))
      SLL_ASSERT(size(x2, 1) == size(y, 1))
      SLL_ASSERT(size(x2, 2) == size(y, 2))

      do i2 = 1, size(x2, 2)
         do i1 = 1, size(x1, 1)
            y(i1, i2) = f_spline_2d__eval_deriv_x1x2(self, x1(i1, i2), x2(i1, i2))
         end do
      end do

   end subroutine s_spline_2d__eval_array_deriv_x1x2

end module sll_m_spline_2d
