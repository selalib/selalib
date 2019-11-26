!> @ingroup splines
!> @brief   Module for 1D splines, linear combination of B-spline functions.
!>
!> @details Here we define a 1D spline type as an element of the linear space
!>          given by the span of the given B-splines (basis functions).
!>          Therefore, initialization of a 1D spline object requires an existing
!>          B-splines object, to which a private (polymorphic) pointer is
!>          associated.
!>          The B-spline coefficients are stored in a public allocatable array;
!>          at initialization the array is allocated to the proper size and all
!>          values are set to zero.
!>          In most situations the B-spline coefficients are not set directly by
!>          the end user, but are computed by some other object (e.g., a Poisson
!>          solver or a spline interpolator).
!>          Various public methods allow the user to evaluate the 1D spline S(x)
!>          and its derivative ∂S(x)/∂x any position x.
!>
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64
  use sll_m_bsplines_base    , only: sll_c_bsplines

  implicit none

  public :: &
    sll_t_spline_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> 1D spline
  type :: sll_t_spline_1d

    real(wp)             , allocatable      :: bcoef(:)
    class(sll_c_bsplines), pointer, private :: bspl => null()

  contains

    procedure :: init             => s_spline_1d__init
    procedure :: free             => s_spline_1d__free
    procedure :: belongs_to_space => f_spline_1d__belongs_to_space
    procedure :: eval             => f_spline_1d__eval
    procedure :: eval_deriv       => f_spline_1d__eval_deriv
    procedure :: eval_array       => s_spline_1d__eval_array
    procedure :: eval_array_deriv => s_spline_1d__eval_array_deriv

  end type sll_t_spline_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief      Initialize 1D spline object as element of span(B-splines)
  !> @param[out] self      1D spline: new element of 1D spline space
  !> @param[in]  bsplines  B-splines: given basis of 1D spline space
  !-----------------------------------------------------------------------------
  subroutine s_spline_1d__init( self, bsplines )

    class(sll_t_spline_1d), intent(  out)         :: self
    class(sll_c_bsplines ), intent(in   ), target :: bsplines

    ! Store pointer to B-splines
    self%bspl => bsplines

    ! Allocate array of spline coefficients: in case of periodic BCs, the last
    ! p coefficients are a periodic copy of the first p ones
    associate( n => bsplines % ncells, p => bsplines % degree )

      allocate( self%bcoef(1:n+p) )

    end associate

    ! Set all coefficients to zero
    self%bcoef = 0.0_f64

  end subroutine s_spline_1d__init

  !-----------------------------------------------------------------------------
  !> @brief        Destroy 1D spline (re-initialization is possible afterwards)
  !> @param[inout] self  1D spline
  !-----------------------------------------------------------------------------
  subroutine s_spline_1d__free( self )

    class(sll_t_spline_1d), intent(inout) :: self

    deallocate( self % bcoef )
    nullify   ( self % bspl  )

  end subroutine s_spline_1d__free

  !-----------------------------------------------------------------------------
  !> @brief     Check if 1D spline belongs to span of given B-splines
  !> @param[in] self      1D spline
  !> @param[in] bsplines  B-splines
  !> @returns             true/false answer
  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d__belongs_to_space( self, bsplines ) result( in_space )

    class(sll_t_spline_1d), intent(in)         :: self
    class(sll_c_bsplines ), intent(in), target :: bsplines
    logical :: in_space

    in_space = associated( self%bspl, bsplines )

  end function f_spline_1d__belongs_to_space

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate value of 1D spline at location x: y=S(x)
  !> @param[in] self  1D spline
  !> @param[in] x     evaluation point
  !> @returns         spline value y=S(x)
  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d__eval( self, x ) result( y )

    class(sll_t_spline_1d), intent(in) :: self
    real(wp)              , intent(in) :: x
    real(wp) :: y

    real(wp) :: values(1:self%bspl%degree+1)
    integer  :: jmin, jmax

    call self % bspl % eval_basis( x, values, jmin )

    jmax = jmin + self%bspl%degree

    y = dot_product( self%bcoef(jmin:jmax), values )

  end function f_spline_1d__eval

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate derivative of 1D spline at location x: y=S'(x)
  !> @param[in] self  1D spline
  !> @param[in] x     evaluation point
  !> @returns         spline derivative y=S'(x)
  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d__eval_deriv( self, x ) result( y )

    class(sll_t_spline_1d), intent(in) :: self
    real(wp)              , intent(in) :: x
    real(wp) :: y

    real(wp) :: derivs(1:self%bspl%degree+1)
    integer  :: jmin, jmax

    call self % bspl % eval_deriv( x, derivs, jmin )

    jmax = jmin + self%bspl%degree

    y = dot_product( self%bcoef(jmin:jmax), derivs )

  end function f_spline_1d__eval_deriv

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate value of 1D spline at all locations in array x
  !> @param[in]  self  1D spline
  !> @param[in]  x     array of evaluation points x[i]
  !> @param[out] y     array of spline values y[i]=S(x[i])
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_1d__eval_array( self, x, y )

    class(sll_t_spline_1d), intent(in   ) :: self
    real(wp)              , intent(in   ) :: x(:)
    real(wp)              , intent(  out) :: y(:)

    integer :: i

    SLL_ASSERT( size(x) == size(y) )

    do i = 1, size(x)
      y(i) = f_spline_1d__eval( self, x(i) )
    end do

  end subroutine s_spline_1d__eval_array

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate derivative of 1D spline at all locations in array x
  !> @param[in]  self  1D spline
  !> @param[in]  x     array of evaluation points x[i]
  !> @param[out] y     array of spline derivatives y[i]=S'(x[i])
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_1d__eval_array_deriv( self, x, y )

    class(sll_t_spline_1d), intent(in   ) :: self
    real(wp)              , intent(in   ) :: x(:)
    real(wp)              , intent(  out) :: y(:)

    integer :: i

    SLL_ASSERT( size(x) == size(y) )

    do i = 1, size(x)
      y(i) = f_spline_1d__eval_deriv( self, x(i) )
    end do

  end subroutine s_spline_1d__eval_array_deriv

end module sll_m_spline_1d
