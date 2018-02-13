module sll_m_polar_mapping_analytical_target
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  implicit none

  public :: sll_t_polar_mapping_analytical_target

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Concrete type, analytical polar mapping
  type, extends(sll_c_polar_mapping_analytical) :: sll_t_polar_mapping_analytical_target

    ! Default values (circle), can be overwritten from 'init' method
    real(wp) :: x0(2) = [ 0.0_wp, 0.0_wp ]
    real(wp) :: d0 = 0.0_wp
    real(wp) :: e0 = 0.0_wp

  contains

    procedure :: init      => s_polar_mapping_analytical_target__init
    procedure :: eval      => f_polar_mapping_analytical_target__eval
    procedure :: jmat      => f_polar_mapping_analytical_target__jmat ! Jacobian matrix
    procedure :: jmat_comp => f_polar_mapping_analytical_target__jmat_comp

  end type sll_t_polar_mapping_analytical_target

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_mapping_analytical_target__init( self, x0, d0, e0 )
    class(sll_t_polar_mapping_analytical_target), intent(inout) :: self
    real(wp), optional                          , intent(in   ) :: x0(2)
    real(wp), optional                          , intent(in   ) :: d0
    real(wp), optional                          , intent(in   ) :: e0

    ! Overwrite parameters
    if ( present( x0 ) ) self % x0 = x0
    if ( present( d0 ) ) self % d0 = d0
    if ( present( e0 ) ) self % e0 = e0

  end subroutine s_polar_mapping_analytical_target__init

  !-----------------------------------------------------------------------------
  SLL_PURE function f_polar_mapping_analytical_target__eval( self, eta ) result( x )
    class(sll_t_polar_mapping_analytical_target), intent(in) :: self
    real(wp)                                    , intent(in) :: eta(2)
    real(wp) :: x(2)

    associate( s => eta(1), t => eta(2), x0 => self % x0, d0 => self % d0, e0 => self % e0 )

      x(1) = x0(1) + (1.0_wp-e0)*s*cos(t) - d0*s**2
      x(2) = x0(2) + (1.0_wp+e0)*s*sin(t)

    end associate

  end function f_polar_mapping_analytical_target__eval

  !-----------------------------------------------------------------------------
  SLL_PURE function f_polar_mapping_analytical_target__jmat( self, eta ) result( jmat )
    class(sll_t_polar_mapping_analytical_target), intent(in) :: self
    real(wp)                                    , intent(in) :: eta(2)
    real(wp) :: jmat(2,2)

    associate( s => eta(1), t => eta(2), d0 => self % d0, e0 => self % e0 )

      ! J_11 = d(x1)/d(eta1)
      ! J_12 = d(x1)/d(eta2)
      ! J_21 = d(x2)/d(eta1)
      ! J_22 = d(x2)/d(eta2)
      jmat(1,1) = -2.0_wp*d0*s + (1.0_wp-e0)*cos(t)
      jmat(1,2) = -(1.0_wp-e0)*s*sin(t)
      jmat(2,1) =  (1.0_wp+e0)*sin(t)
      jmat(2,2) =  (1.0_wp+e0)*s*cos(t)

    end associate

  end function f_polar_mapping_analytical_target__jmat

  !-----------------------------------------------------------------------------
  SLL_PURE function f_polar_mapping_analytical_target__jmat_comp( self, eta ) result( jmat_comp )
    class(sll_t_polar_mapping_analytical_target), intent(in) :: self
    real(wp)                                    , intent(in) :: eta(2)
    real(wp) :: jmat_comp(2,2)

    associate( s => eta(1), t => eta(2), d0 => self % d0, e0 => self % e0 )

      jmat_comp(1,1) = (1.0_wp+e0)
      jmat_comp(1,2) = 2.0_wp*d0*s*sin(t)
      jmat_comp(2,1) = 0.0_wp
      jmat_comp(2,2) = (1.0_wp-e0)-2.0_wp*d0*s*cos(t)

      jmat_comp = jmat_comp / ( (1.0_wp+e0)*((1.0_wp-e0)-2.0_wp*d0*s*cos(t)) )

    end associate

  end function f_polar_mapping_analytical_target__jmat_comp

end module sll_m_polar_mapping_analytical_target
