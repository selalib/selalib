module sll_m_polar_mapping_analytical_czarny
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  implicit none

  public :: sll_t_polar_mapping_analytical_czarny

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Concrete type, analytical polar mapping
  type, extends(sll_c_polar_mapping_analytical) :: sll_t_polar_mapping_analytical_czarny

    ! Default values (circle), can be overwritten from 'init' method
    real(wp) :: x0(2) = [ 0.0_wp, 0.0_wp ]
    real(wp) :: b = 1.0_wp
    real(wp) :: e = 0.0_wp

  contains

    procedure :: init => s_polar_mapping_analytical_czarny__init
    procedure :: eval => f_polar_mapping_analytical_czarny__eval
    procedure :: jmat => f_polar_mapping_analytical_czarny__jmat ! Jacobian matrix

  end type sll_t_polar_mapping_analytical_czarny

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_mapping_analytical_czarny__init( self, x0, b, e )
    class(sll_t_polar_mapping_analytical_czarny), intent(inout) :: self
    real(wp), optional                          , intent(in   ) :: x0(2)
    real(wp), optional                          , intent(in   ) :: b
    real(wp), optional                          , intent(in   ) :: e

    ! Overwrite parameters
    if ( present( x0 ) ) self % x0 = x0
    if ( present( b  ) ) self % b  = b
    if ( present( e  ) ) self % e  = e

  end subroutine s_polar_mapping_analytical_czarny__init

  !-----------------------------------------------------------------------------
  SLL_PURE function f_polar_mapping_analytical_czarny__eval( self, eta ) result( x )
    class(sll_t_polar_mapping_analytical_czarny), intent(in) :: self
    real(wp)                                    , intent(in) :: eta(2)
    real(wp) :: x(2)

    associate( s => eta(1), t => eta(2), x0 => self % x0, b => self % b, e => self % e )

      ! TODO: check if x0(1) can be added to x(1) as well
      x(1) = (1.0_wp-sqrt( 1.0_wp + e*(e+2.0_wp*s*cos(t)) ))/e
      x(2) = x0(2) + b*s*sin(t)/((1.0_wp+e*x(1))*sqrt( 1.0_wp-e**2*0.25_wp ))

    end associate

  end function f_polar_mapping_analytical_czarny__eval

  !-----------------------------------------------------------------------------
  SLL_PURE function f_polar_mapping_analytical_czarny__jmat( self, eta ) result( jmat )
    class(sll_t_polar_mapping_analytical_czarny), intent(in) :: self
    real(wp)                                    , intent(in) :: eta(2)
    real(wp) :: jmat(2,2)

    real(wp) :: tmp1, tmp2

    associate( s => eta(1), t => eta(2), b => self % b, e => self % e )

      tmp1 = sqrt( 1.0_wp + e*(e+2.0_wp*s*cos(t)) )
      tmp2 = sqrt( 1.0_wp-e**2*0.25_wp )

      ! J_11 = d(x1)/d(eta1)
      ! J_12 = d(x1)/d(eta2)
      ! J_21 = d(x2)/d(eta1)
      ! J_22 = d(x2)/d(eta2)
      jmat(1,1) = -cos(t)/tmp1
      jmat(1,2) = s*sin(t)/tmp1
      jmat(2,1) = b*sin(t)/((2.0_wp-tmp1)*tmp2) + e*b*s*sin(t)*cos(t)/(tmp1*tmp2*(2.0_wp-tmp1)**2)
      jmat(2,2) = b*s*cos(t)/((2.0_wp-tmp1)*tmp2) - e*b*s**2*sin(t)**2/(tmp1*tmp2*(2.0_wp-tmp1)**2)

    end associate

  end function f_polar_mapping_analytical_czarny__jmat

end module sll_m_polar_mapping_analytical_czarny
