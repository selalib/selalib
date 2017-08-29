module sll_m_spline_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64
  use sll_m_bsplines_base    , only: sll_c_bsplines

  implicit none

  public :: &
    sll_t_spline_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> 1D spline
  type :: sll_t_spline_2d

    real(wp)             , allocatable      :: bcoef(:,:)
    class(sll_c_bsplines), pointer, private :: bspl1 => null()
    class(sll_c_bsplines), pointer, private :: bspl2 => null()

  contains

    procedure :: init                => s_spline_2d__init
    procedure :: free                => s_spline_2d__free
    procedure :: eval                => f_spline_2d__eval
    procedure :: eval_deriv_x1       => f_spline_2d__eval_deriv_x1
    procedure :: eval_deriv_x2       => f_spline_2d__eval_deriv_x2
    procedure :: eval_array          => s_spline_2d__eval_array
    procedure :: eval_array_deriv_x1 => s_spline_2d__eval_array_deriv_x1
    procedure :: eval_array_deriv_x2 => s_spline_2d__eval_array_deriv_x2

  end type sll_t_spline_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_spline_2d__init( self, bsplines_x1, bsplines_x2 )

    class(sll_t_spline_2d), intent(  out)         :: self
    class(sll_c_bsplines) , intent(in   ), target :: bsplines_x1
    class(sll_c_bsplines) , intent(in   ), target :: bsplines_x2

    ! Store pointer to B-splines
    self%bspl1 => bsplines_x1
    self%bspl2 => bsplines_x2

    ! Allocate array of spline coefficients
    ! in case of periodic BCs, a larger array of coefficients is used in order
    ! to avoid a loop with calls to the "mod( , )" function at evaluation.
    associate( n1 => bsplines_x1 % nbasis, &
               n2 => bsplines_x2 % nbasis, &
               g1 => merge( 1 + bsplines_x1 % degree/2, 0, bsplines_x1 % periodic ), &
               g2 => merge( 1 + bsplines_x2 % degree/2, 0, bsplines_x2 % periodic ) )

      allocate( self%bcoef(1-g1:n1+g1,1-g2:n2+g2) )
      
    end associate

    ! Set all coefficients to zero
    self%bcoef = 0.0_f64

  end subroutine s_spline_2d__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_2d__free( self )

    class(sll_t_spline_2d), intent(inout) :: self

    deallocate( self % bcoef )
    nullify   ( self % bspl1 )
    nullify   ( self % bspl2 )

  end subroutine s_spline_2d__free

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_2d__eval( self, x1, x2 ) result( y )

    class(sll_t_spline_2d), intent(in) :: self
    real(wp)              , intent(in) :: x1
    real(wp)              , intent(in) :: x2
    real(wp) :: y

    integer :: jmin(2)
    integer :: jmax(2)
    integer :: k1, k2

    ! Automatic arrays
    real(wp) :: values1(1:self%bspl1%degree+1)
    real(wp) :: values2(1:self%bspl2%degree+1)

    ! Compute arrays v1 and v2 of B-spline values
    call self % bspl1 % eval_basis( x1, values1, jmin(1) )
    call self % bspl2 % eval_basis( x2, values2, jmin(2) )

    jmax(1) = jmin(1) + self%bspl1%degree
    jmax(2) = jmin(2) + self%bspl2%degree

    ! Determine (matrix) block C of B-spline coefficients to be used
    ! and compute scalar product <v1,v2> = (v1^T)*C*v2
    associate( bcoef => self%bcoef (jmin(1):jmax(1), jmin(2):jmax(2)) )
      y = 0.0_f64
      do k2 = 1, 1+self%bspl2%degree
        do k1 = 1, 1+self%bspl1%degree
          y  = y + bcoef(k1,k2) * values1(k1) * values2(k2)
        end do
      end do
    end associate

  end function f_spline_2d__eval

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_2d__eval_deriv_x1( self, x1, x2 ) result( y )

    class(sll_t_spline_2d), intent(in) :: self
    real(wp)              , intent(in) :: x1
    real(wp)              , intent(in) :: x2
    real(wp) :: y

    integer :: jmin(2)
    integer :: jmax(2)
    integer :: k1, k2

    ! Automatic arrays
    real(wp) :: derivs1(1:self%bspl1%degree+1)
    real(wp) :: values2(1:self%bspl2%degree+1)

    ! Compute arrays d1 and v2 of B-spline derivatives/values
    call self % bspl1 % eval_deriv( x1, derivs1, jmin(1) )
    call self % bspl2 % eval_basis( x2, values2, jmin(2) )

    jmax(1) = jmin(1) + self%bspl1%degree
    jmax(2) = jmin(2) + self%bspl2%degree

    ! Determine (matrix) block C of B-spline coefficients to be used
    ! and compute scalar product <d1,v2> = (d1^T)*C*v2
    associate( bcoef => self%bcoef (jmin(1):jmax(1), jmin(2):jmax(2)) )
      y = 0.0_f64
      do k2 = 1, 1+self%bspl2%degree
        do k1 = 1, 1+self%bspl1%degree
          y  = y + bcoef(k1,k2) * derivs1(k1) * values2(k2)
        end do
      end do
    end associate

  end function f_spline_2d__eval_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_2d__eval_deriv_x2( self, x1, x2 ) result( y )

    class(sll_t_spline_2d), intent(in) :: self
    real(wp)              , intent(in) :: x1
    real(wp)              , intent(in) :: x2
    real(wp) :: y

    integer :: jmin(2)
    integer :: jmax(2)
    integer :: k1, k2

    ! Automatic arrays
    real(wp) :: values1(1:self%bspl1%degree+1)
    real(wp) :: derivs2(1:self%bspl2%degree+1)

    ! Compute arrays v1 and d2 of B-spline values/derivatives
    call self % bspl1 % eval_basis( x1, values1, jmin(1) )
    call self % bspl2 % eval_deriv( x2, derivs2, jmin(2) )

    jmax(1) = jmin(1) + self%bspl1%degree
    jmax(2) = jmin(2) + self%bspl2%degree

    ! Determine (matrix) block C of B-spline coefficients to be used
    ! and compute scalar product <v1,d2> = (v1^T)*C*d2
    associate( bcoef => self%bcoef (jmin(1):jmax(1), jmin(2):jmax(2)) )
      y = 0.0_f64
      do k2 = 1, 1+self%bspl2%degree
        do k1 = 1, 1+self%bspl1%degree
          y  = y + bcoef(k1,k2) * values1(k1) * derivs2(k2)
        end do
      end do
    end associate

  end function f_spline_2d__eval_deriv_x2

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_2d__eval_array( self, x1, x2, y )

    class(sll_t_spline_2d), intent(in   ) :: self
    real(wp)              , intent(in   ) :: x1(:,:)
    real(wp)              , intent(in   ) :: x2(:,:)
    real(wp)              , intent(  out) :: y (:,:)

    integer :: i1, i2

    SLL_ASSERT( size(x1,1) == size(y,1) )
    SLL_ASSERT( size(x1,2) == size(y,2) )
    SLL_ASSERT( size(x2,1) == size(y,1) )
    SLL_ASSERT( size(x2,2) == size(y,2) )

    do i2 = 1, size(x2,2)
      do i1 = 1, size(x1,1)
        y(i1,i2) = f_spline_2d__eval( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine s_spline_2d__eval_array

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_2d__eval_array_deriv_x1( self, x1, x2, y )

    class(sll_t_spline_2d), intent(in   ) :: self
    real(wp)              , intent(in   ) :: x1(:,:)
    real(wp)              , intent(in   ) :: x2(:,:)
    real(wp)              , intent(  out) :: y (:,:)

    integer :: i1, i2

    SLL_ASSERT( size(x1,1) == size(y,1) )
    SLL_ASSERT( size(x1,2) == size(y,2) )
    SLL_ASSERT( size(x2,1) == size(y,1) )
    SLL_ASSERT( size(x2,2) == size(y,2) )

    do i2 = 1, size(x2,2)
      do i1 = 1, size(x1,1)
        y(i1,i2) = f_spline_2d__eval_deriv_x1( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine s_spline_2d__eval_array_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_2d__eval_array_deriv_x2( self, x1, x2, y )

    class(sll_t_spline_2d), intent(in   ) :: self
    real(wp)              , intent(in   ) :: x1(:,:)
    real(wp)              , intent(in   ) :: x2(:,:)
    real(wp)              , intent(  out) :: y (:,:)

    integer :: i1, i2

    SLL_ASSERT( size(x1,1) == size(y,1) )
    SLL_ASSERT( size(x1,2) == size(y,2) )
    SLL_ASSERT( size(x2,1) == size(y,1) )
    SLL_ASSERT( size(x2,2) == size(y,2) )

    do i2 = 1, size(x2,2)
      do i1 = 1, size(x1,1)
        y(i1,i2) = f_spline_2d__eval_deriv_x2( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine s_spline_2d__eval_array_deriv_x2

end module sll_m_spline_2d_new
