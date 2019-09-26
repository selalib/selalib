module sll_m_polar_bsplines_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_bsplines_base, only: sll_c_bsplines

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

  implicit none

  public :: sll_t_polar_bsplines_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Type containing new 2D polar basis functions
  type :: sll_t_polar_bsplines_2d

    ! Pointer to discrete IGA mapping (private member)
    type(sll_t_singular_mapping_discrete), pointer, private :: mapping => null()

    ! Pole of the singularity
    real(wp) :: pole(2)

    ! Parameter to compute triangle vertices and barycentric coordinates
    real(wp) :: tau

    ! Triangle vertices
    real(wp) :: T0(2)
    real(wp) :: T1(2)
    real(wp) :: T2(2)

    ! New basis functions
    type(sll_t_spline_2d) :: N0
    type(sll_t_spline_2d) :: N1
    type(sll_t_spline_2d) :: N2

  contains

    procedure :: init    => s_polar_bsplines_2d__init
    procedure :: eval_l0 => f_polar_bsplines_2d__eval_l0
    procedure :: eval_l1 => f_polar_bsplines_2d__eval_l1
    procedure :: eval_l2 => f_polar_bsplines_2d__eval_l2
    procedure :: free    => s_polar_bsplines_2d__free

  end type sll_t_polar_bsplines_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_bsplines_2d__init( self, spline_basis_eta1, spline_basis_eta2, mapping )
    class(sll_t_polar_bsplines_2d)               , intent(inout) :: self
    class(sll_c_bsplines)                        , intent(in   ) :: spline_basis_eta1
    class(sll_c_bsplines)                        , intent(in   ) :: spline_basis_eta2
    type(sll_t_singular_mapping_discrete), target, intent(in   ) :: mapping ! mapping must be discrete

    real(wp) :: tau1, tau2, tau3
    integer  :: nb1, nb2

    ! Assign pointer to discrete IGA mapping
    self % mapping => mapping

    ! Store pole
    self % pole = mapping % pole()

    ! Store locally number of basis functions
    nb1 = spline_basis_eta1 % nbasis
    nb2 = spline_basis_eta2 % nbasis

    ! Compute tau from mapping control points
    associate( x0   => self % pole(1), &
               y0   => self % pole(2), &
               c_x1 => mapping % spline_2d_x1 % bcoef(2,1:nb2), &
               c_x2 => mapping % spline_2d_x2 % bcoef(2,1:nb2) )

      tau1 = maxval( - 2.0_wp * (c_x1-x0) )
      tau2 = maxval( c_x1 - x0 - sqrt(3.0_wp) * (c_x2-y0) )
      tau3 = maxval( c_x1 - x0 + sqrt(3.0_wp) * (c_x2-y0) )

      self % tau = max( tau1, tau2, tau3 )

      ! Compute triangle vertices
      self % T0(1) = x0 + self % tau
      self % T0(2) = y0
      self % T1(1) = x0 - self % tau * 0.5_wp
      self % T1(2) = y0 + self % tau * 0.5_wp * sqrt(3.0_wp)
      self % T2(1) = x0 - self % tau * 0.5_wp
      self % T2(2) = y0 - self % tau * 0.5_wp * sqrt(3.0_wp)

    end associate

    ! Initialize new basis functions as 2D splines
    call self % N0 % init( spline_basis_eta1, spline_basis_eta2 )
    call self % N1 % init( spline_basis_eta1, spline_basis_eta2 )
    call self % N2 % init( spline_basis_eta1, spline_basis_eta2 )

    ! Set to zero all coefficients
    self % N0 % bcoef(1:nb1,1:nb2) = 0.0_wp
    self % N1 % bcoef(1:nb1,1:nb2) = 0.0_wp
    self % N2 % bcoef(1:nb1,1:nb2) = 0.0_wp

    ! Set coefficients e_0j (i=1)
    self % N0 % bcoef(1,1:nb2) = 1.0_wp / 3.0_wp
    self % N1 % bcoef(1,1:nb2) = 1.0_wp / 3.0_wp
    self % N2 % bcoef(1,1:nb2) = 1.0_wp / 3.0_wp

    ! Set coefficients e_1j (i=2)
    associate( c_x1 => mapping % spline_2d_x1 % bcoef(2,1:nb2), &
               c_x2 => mapping % spline_2d_x2 % bcoef(2,1:nb2) )

      self % N0 % bcoef(2,1:nb2) = self % eval_l0( c_x1, c_x2 )
      self % N1 % bcoef(2,1:nb2) = self % eval_l1( c_x1, c_x2 )
      self % N2 % bcoef(2,1:nb2) = self % eval_l2( c_x1, c_x2 )

    end associate

  end subroutine s_polar_bsplines_2d__init

  !-----------------------------------------------------------------------------
  SLL_PURE elemental function f_polar_bsplines_2d__eval_l0( self, x, y ) result( l0 )
    class(sll_t_polar_bsplines_2d), intent(in) :: self
    real(wp)                      , intent(in) :: x
    real(wp)                      , intent(in) :: y 
    real(wp) :: l0

    associate( x0 => self % pole(1), y0 => self % pole(2), tau => self % tau )

      l0 = 1.0_wp / 3.0_wp + 2.0_wp * (x-x0) / ( 3.0_wp * tau )

    end associate

  end function f_polar_bsplines_2d__eval_l0

  !-----------------------------------------------------------------------------
  SLL_PURE elemental function f_polar_bsplines_2d__eval_l1( self, x, y ) result( l1 )
    class(sll_t_polar_bsplines_2d), intent(in) :: self
    real(wp)                      , intent(in) :: x
    real(wp)                      , intent(in) :: y
    real(wp) :: l1

    associate( x0 => self % pole(1), y0 => self % pole(2), tau => self % tau )

      l1 = 1.0_wp / 3.0_wp - (x-x0) / ( 3.0_wp * tau ) + (y-y0) / ( sqrt(3.0_wp) * tau )

    end associate

  end function f_polar_bsplines_2d__eval_l1

  !-----------------------------------------------------------------------------
  SLL_PURE elemental function f_polar_bsplines_2d__eval_l2( self, x, y ) result( l2 )
    class(sll_t_polar_bsplines_2d), intent(in) :: self
    real(wp)                      , intent(in) :: x
    real(wp)                      , intent(in) :: y
    real(wp) :: l2

    associate( x0 => self % pole(1), y0 => self % pole(2), tau => self % tau )

      l2 = 1.0_wp / 3.0_wp - (x-x0) / ( 3.0_wp * tau ) - (y-y0) / ( sqrt(3.0_wp) * tau )

    end associate

  end function f_polar_bsplines_2d__eval_l2

  !-----------------------------------------------------------------------------
  subroutine s_polar_bsplines_2d__free( self )
    class(sll_t_polar_bsplines_2d), intent(inout) :: self

    call self % N0 % free()
    call self % N1 % free()
    call self % N2 % free()

  end subroutine s_polar_bsplines_2d__free

end module sll_m_polar_bsplines_2d
