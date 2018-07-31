module sll_m_jacobian_2d_pseudo_cartesian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_pi

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  implicit none

  public :: sll_t_jacobian_2d_pseudo_cartesian

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_jacobian_2d_pseudo_cartesian

    ! Can be overwritten at initialization
    real(wp) :: eps = 1.0e-12_wp

    type(sll_t_polar_mapping_iga), pointer :: mapping => null()

  contains

    procedure :: init => s_jacobian_2d_pseudo_cartesian__init
    procedure :: eval => s_jacobian_2d_pseudo_cartesian__eval
    procedure :: free => s_jacobian_2d_pseudo_cartesian__free

  end type sll_t_jacobian_2d_pseudo_cartesian

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_jacobian_2d_pseudo_cartesian__init( self, mapping, eps )
    class(sll_t_jacobian_2d_pseudo_cartesian), intent(inout) :: self
    type (sll_t_polar_mapping_iga), target   , intent(in   ) :: mapping
    real(wp), optional                       , intent(in   ) :: eps

    self % mapping => mapping

    if ( present( eps ) ) self % eps = eps

  end subroutine s_jacobian_2d_pseudo_cartesian__init

  !-----------------------------------------------------------------------------
  SLL_PURE function s_jacobian_2d_pseudo_cartesian__eval( self, eta ) result( jmat )
    class(sll_t_jacobian_2d_pseudo_cartesian), intent(in) :: self
    real(wp)                                 , intent(in) :: eta(2)
    real(wp) :: jmat(2,2) ! inverse of Jacobian

    real(wp) :: eps, jmat_pole(2,2), jmat_eps(2,2)

    eps = self % eps

    if ( eta(1) == 0.0_wp ) then

      jmat(2,2) =   self % mapping % spline_2d_x1 % eval_deriv_x1  ( eta(1), eta(2) ) * cos( eta(2) ) &
                  - self % mapping % spline_2d_x1 % eval_deriv_x1x2( eta(1), eta(2) ) * sin( eta(2) )

      jmat(1,2) = - self % mapping % spline_2d_x1 % eval_deriv_x1  ( eta(1), eta(2) ) * sin( eta(2) ) &
                  - self % mapping % spline_2d_x1 % eval_deriv_x1x2( eta(1), eta(2) ) * cos( eta(2) )

      jmat(2,1) = - self % mapping % spline_2d_x2 % eval_deriv_x1  ( eta(1), eta(2) ) * cos( eta(2) ) &
                  + self % mapping % spline_2d_x2 % eval_deriv_x1x2( eta(1), eta(2) ) * sin( eta(2) )

      jmat(1,1) =   self % mapping % spline_2d_x2 % eval_deriv_x1  ( eta(1), eta(2) ) * sin( eta(2) ) &
                  + self % mapping % spline_2d_x2 % eval_deriv_x1x2( eta(1), eta(2) ) * cos( eta(2) )

      jmat(:,:) = jmat(:,:) / ( jmat(1,1)*jmat(2,2) - jmat(1,2)*jmat(2,1) )

    else if ( 0.0_wp < eta(1) .and. eta(1) < eps ) then

      ! Jacobian at s = 0

      jmat_pole(2,2) =   self % mapping % spline_2d_x1 % eval_deriv_x1  ( eta(1), eta(2) ) * cos( eta(2) ) &
                       - self % mapping % spline_2d_x1 % eval_deriv_x1x2( eta(1), eta(2) ) * sin( eta(2) )

      jmat_pole(1,2) = - self % mapping % spline_2d_x1 % eval_deriv_x1  ( eta(1), eta(2) ) * sin( eta(2) ) &
                       - self % mapping % spline_2d_x1 % eval_deriv_x1x2( eta(1), eta(2) ) * cos( eta(2) )

      jmat_pole(2,1) = - self % mapping % spline_2d_x2 % eval_deriv_x1  ( eta(1), eta(2) ) * cos( eta(2) ) &
                       + self % mapping % spline_2d_x2 % eval_deriv_x1x2( eta(1), eta(2) ) * sin( eta(2) )

      jmat_pole(1,1) =   self % mapping % spline_2d_x2 % eval_deriv_x1  ( eta(1), eta(2) ) * sin( eta(2) ) &
                       + self % mapping % spline_2d_x2 % eval_deriv_x1x2( eta(1), eta(2) ) * cos( eta(2) )

      jmat_pole(:,:) = jmat_pole(:,:) / ( jmat_pole(1,1)*jmat_pole(2,2) - jmat_pole(1,2)*jmat_pole(2,1) )

      ! Jacobian at s = eps

      jmat_eps(2,2) =   self % mapping % spline_2d_x1 % eval_deriv_x1( eps, eta(2) ) * cos( eta(2) ) &
                      - self % mapping % spline_2d_x1 % eval_deriv_x2( eps, eta(2) ) * sin( eta(2) ) / eps

      jmat_eps(1,2) = - self % mapping % spline_2d_x1 % eval_deriv_x1( eps, eta(2) ) * sin( eta(2) ) &
                      - self % mapping % spline_2d_x1 % eval_deriv_x2( eps, eta(2) ) * cos( eta(2) ) / eps

      jmat_eps(2,1) = - self % mapping % spline_2d_x2 % eval_deriv_x1( eps, eta(2) ) * cos( eta(2) ) &
                      + self % mapping % spline_2d_x2 % eval_deriv_x2( eps, eta(2) ) * sin( eta(2) ) / eps

      jmat_eps(1,1) =   self % mapping % spline_2d_x2 % eval_deriv_x1( eps, eta(2) ) * sin( eta(2) ) &
                      + self % mapping % spline_2d_x2 % eval_deriv_x2( eps, eta(2) ) * cos( eta(2) ) / eps

      jmat_eps(:,:) = jmat_eps(:,:) / ( jmat_eps(1,1)*jmat_eps(2,2) - jmat_eps(1,2)*jmat_eps(2,1) )

      ! Linear interpolation between s = 0 and s = eps

      jmat(1,1) = (1.0_wp-eta(1)/eps)*jmat_pole(1,1) + eta(1)/eps*jmat_eps(1,1)
      jmat(1,2) = (1.0_wp-eta(1)/eps)*jmat_pole(1,2) + eta(1)/eps*jmat_eps(1,2)
      jmat(2,1) = (1.0_wp-eta(1)/eps)*jmat_pole(2,1) + eta(1)/eps*jmat_eps(2,1)
      jmat(2,2) = (1.0_wp-eta(1)/eps)*jmat_pole(2,2) + eta(1)/eps*jmat_eps(2,2)

    else if ( eps <= eta(1) ) then

      jmat(2,2) =   self % mapping % spline_2d_x1 % eval_deriv_x1( eta(1), eta(2) ) * cos( eta(2) ) &
                  - self % mapping % spline_2d_x1 % eval_deriv_x2( eta(1), eta(2) ) * sin( eta(2) ) / eta(1)

      jmat(1,2) = - self % mapping % spline_2d_x1 % eval_deriv_x1( eta(1), eta(2) ) * sin( eta(2) ) &
                  - self % mapping % spline_2d_x1 % eval_deriv_x2( eta(1), eta(2) ) * cos( eta(2) ) / eta(1)

      jmat(2,1) = - self % mapping % spline_2d_x2 % eval_deriv_x1( eta(1), eta(2) ) * cos( eta(2) ) &
                  + self % mapping % spline_2d_x2 % eval_deriv_x2( eta(1), eta(2) ) * sin( eta(2) ) / eta(1)

      jmat(1,1) =   self % mapping % spline_2d_x2 % eval_deriv_x1( eta(1), eta(2) ) * sin( eta(2) ) &
                  + self % mapping % spline_2d_x2 % eval_deriv_x2( eta(1), eta(2) ) * cos( eta(2) ) / eta(1)

      jmat(:,:) = jmat(:,:) / ( jmat(1,1)*jmat(2,2) - jmat(1,2)*jmat(2,1) )

    end if

  end function s_jacobian_2d_pseudo_cartesian__eval

  !-----------------------------------------------------------------------------
  subroutine s_jacobian_2d_pseudo_cartesian__free( self )
    class(sll_t_jacobian_2d_pseudo_cartesian), intent(inout) :: self

    nullify( self % mapping )

  end subroutine s_jacobian_2d_pseudo_cartesian__free

end module sll_m_jacobian_2d_pseudo_cartesian
