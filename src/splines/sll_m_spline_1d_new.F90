module sll_m_spline_1d_new
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

    procedure :: init                => s_spline_1d__init
    procedure :: free                => s_spline_1d__free
    procedure :: eval                => f_spline_1d__eval
    procedure :: eval_deriv          => f_spline_1d__eval_deriv
    procedure :: eval_array          => s_spline_1d__eval_array
    procedure :: eval_array_deriv    => s_spline_1d__eval_array_deriv

  end type sll_t_spline_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_spline_1d__init( self, bsplines )

    class(sll_t_spline_1d), intent(  out)         :: self
    class(sll_c_bsplines) , intent(in   ), target :: bsplines

    ! Store pointer to B-splines
    self%bspl => bsplines

    ! Allocate array of spline coefficients
    ! in case of periodic BCs, a larger array of coefficients is used in order
    ! to avoid a loop with calls to the "mod( , )" function at evaluation.
    associate( n => bsplines % nbasis, &
               g => merge( 1 + bsplines % degree/2, 0, bsplines % periodic ) )

      allocate( self%bcoef(1-g:n+g) )

    end associate

    ! Set all coefficients to zero
    self%bcoef = 0.0_f64

  end subroutine s_spline_1d__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_1d__free( self )

    class(sll_t_spline_1d), intent(inout) :: self

    deallocate( self % bcoef )
    nullify   ( self % bspl  )

  end subroutine s_spline_1d__free

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

end module sll_m_spline_1d_new
