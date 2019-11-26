! This module implements the analytical mapping given by
! equation (2) of https://doi.org/10.1016/j.jcp.2019.108889:
!
! x(s,theta)=x0+(1-kappa)*s*cos(theta)-Delta*s**2
! y(s,theta)=y0+(1+kappa)*s*sin(theta)

module sll_m_singular_mapping_analytic_target
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_singular_mapping_analytic, only: sll_c_singular_mapping_analytic

  implicit none

  public :: sll_t_singular_mapping_analytic_target

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Concrete type, analytical singular mapping
  type, extends(sll_c_singular_mapping_analytic) :: sll_t_singular_mapping_analytic_target

    ! Default values (circle), can be overwritten from 'init' method
    real(wp) :: x0(2) = [ 0.0_wp, 0.0_wp ]
    real(wp) :: Delta = 0.0_wp
    real(wp) :: kappa = 0.0_wp

  contains

    procedure :: init      => s_singular_mapping_analytic_target__init
    procedure :: eval      => f_singular_mapping_analytic_target__eval
    procedure :: jmat      => f_singular_mapping_analytic_target__jmat ! Jacobian matrix
    procedure :: jmat_comp => f_singular_mapping_analytic_target__jmat_comp

  end type sll_t_singular_mapping_analytic_target

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_singular_mapping_analytic_target__init( self, x0, Delta, kappa )
    class(sll_t_singular_mapping_analytic_target), intent(inout) :: self
    real(wp), optional                           , intent(in   ) :: x0(2)
    real(wp), optional                           , intent(in   ) :: Delta
    real(wp), optional                           , intent(in   ) :: kappa

    ! Overwrite parameters
    if ( present( x0 ) ) self % x0 = x0
    if ( present( Delta ) ) self % Delta = Delta
    if ( present( kappa ) ) self % kappa = kappa

  end subroutine s_singular_mapping_analytic_target__init

  !-----------------------------------------------------------------------------
  SLL_PURE function f_singular_mapping_analytic_target__eval( self, eta ) result( x )
    class(sll_t_singular_mapping_analytic_target), intent(in) :: self
    real(wp)                                     , intent(in) :: eta(2)
    real(wp) :: x(2)

    associate( s => eta(1), t => eta(2), x0 => self % x0, Delta => self % Delta, kappa => self % kappa )

      x(1) = x0(1) + (1.0_wp-kappa)*s*cos(t) - Delta*s**2
      x(2) = x0(2) + (1.0_wp+kappa)*s*sin(t)

    end associate

  end function f_singular_mapping_analytic_target__eval

  !-----------------------------------------------------------------------------
  SLL_PURE function f_singular_mapping_analytic_target__jmat( self, eta ) result( jmat )
    class(sll_t_singular_mapping_analytic_target), intent(in) :: self
    real(wp)                                     , intent(in) :: eta(2)
    real(wp) :: jmat(2,2)

    associate( s => eta(1), t => eta(2), Delta => self % Delta, kappa => self % kappa )

      ! J_11 = d(x1)/d(eta1)
      ! J_12 = d(x1)/d(eta2)
      ! J_21 = d(x2)/d(eta1)
      ! J_22 = d(x2)/d(eta2)
      jmat(1,1) = -2.0_wp*Delta*s + (1.0_wp-kappa)*cos(t)
      jmat(1,2) = -(1.0_wp-kappa)*s*sin(t)
      jmat(2,1) =  (1.0_wp+kappa)*sin(t)
      jmat(2,2) =  (1.0_wp+kappa)*s*cos(t)

    end associate

  end function f_singular_mapping_analytic_target__jmat

  !-----------------------------------------------------------------------------
  SLL_PURE function f_singular_mapping_analytic_target__jmat_comp( self, eta ) result( jmat_comp )
    class(sll_t_singular_mapping_analytic_target), intent(in) :: self
    real(wp)                                     , intent(in) :: eta(2)
    real(wp) :: jmat_comp(2,2)

    associate( s => eta(1), t => eta(2), Delta => self % Delta, kappa => self % kappa )

      jmat_comp(1,1) = (1.0_wp+kappa)
      jmat_comp(1,2) = 2.0_wp*Delta*s*sin(t)
      jmat_comp(2,1) = 0.0_wp
      jmat_comp(2,2) = (1.0_wp-kappa)-2.0_wp*Delta*s*cos(t)

      jmat_comp = jmat_comp / ( (1.0_wp+kappa)*((1.0_wp-kappa)-2.0_wp*Delta*s*cos(t)) )

    end associate

  end function f_singular_mapping_analytic_target__jmat_comp

end module sll_m_singular_mapping_analytic_target
