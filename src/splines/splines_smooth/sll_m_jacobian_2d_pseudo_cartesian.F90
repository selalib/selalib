module sll_m_jacobian_2d_pseudo_cartesian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_pi

  use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

  implicit none

  public :: sll_t_jacobian_2d_pseudo_cartesian

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_jacobian_2d_pseudo_cartesian

    ! Can be overwritten at initialization
    real(wp) :: eps = 1.0e-12_wp

    type(sll_t_singular_mapping_discrete), pointer :: mapping => null()

    real(wp) :: jmat_pole(2,2)

  contains

    procedure :: init => s_jacobian_2d_pseudo_cartesian__init
    procedure :: pole => s_jacobian_2d_pseudo_cartesian__pole
    procedure :: eval => s_jacobian_2d_pseudo_cartesian__eval
    procedure :: free => s_jacobian_2d_pseudo_cartesian__free

  end type sll_t_jacobian_2d_pseudo_cartesian

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_jacobian_2d_pseudo_cartesian__init( self, mapping, eps )
    class(sll_t_jacobian_2d_pseudo_cartesian)     , intent(inout) :: self
    type (sll_t_singular_mapping_discrete), target, intent(in   ) :: mapping
    real(wp), optional                            , intent(in   ) :: eps

    self % mapping => mapping

    if ( present( eps ) ) self % eps = eps

  end subroutine s_jacobian_2d_pseudo_cartesian__init

  !-----------------------------------------------------------------------------
  subroutine s_jacobian_2d_pseudo_cartesian__pole( self, eta2_array )
    class(sll_t_jacobian_2d_pseudo_cartesian), intent(inout) :: self
    real(wp)                                 , intent(in   ) :: eta2_array(:)

    integer  :: i2
    real(wp) :: s, theta
    real(wp) :: jmat_temp(2,2)

    self % jmat_pole(:,:) = 0.0_wp

    s = 0.0_wp

    do i2 = 1, size( eta2_array )

      theta = eta2_array(i2)

      jmat_temp(2,2) =   self % mapping % spline_2d_x1 % eval_deriv_x1  ( s, theta ) * cos( theta ) &
                       - self % mapping % spline_2d_x1 % eval_deriv_x1x2( s, theta ) * sin( theta )

      jmat_temp(1,2) = - self % mapping % spline_2d_x1 % eval_deriv_x1  ( s, theta ) * sin( theta ) &
                       - self % mapping % spline_2d_x1 % eval_deriv_x1x2( s, theta ) * cos( theta )

      jmat_temp(2,1) = - self % mapping % spline_2d_x2 % eval_deriv_x1  ( s, theta ) * cos( theta ) &
                       + self % mapping % spline_2d_x2 % eval_deriv_x1x2( s, theta ) * sin( theta )

      jmat_temp(1,1) =   self % mapping % spline_2d_x2 % eval_deriv_x1  ( s, theta ) * sin( theta ) &
                       + self % mapping % spline_2d_x2 % eval_deriv_x1x2( s, theta ) * cos( theta )

      jmat_temp(:,:) = jmat_temp(:,:) / ( jmat_temp(1,1) * jmat_temp(2,2) - jmat_temp(1,2) * jmat_temp(2,1) )

      self % jmat_pole(:,:) = self % jmat_pole(:,:) + jmat_temp(:,:)

    end do

    ! Take average over all values of theta
    self % jmat_pole(:,:) = self % jmat_pole(:,:) / size( eta2_array )

  end subroutine s_jacobian_2d_pseudo_cartesian__pole

  !-----------------------------------------------------------------------------
  SLL_PURE function s_jacobian_2d_pseudo_cartesian__eval( self, eta ) result( jmat )
    class(sll_t_jacobian_2d_pseudo_cartesian), intent(in) :: self
    real(wp)                                 , intent(in) :: eta(2)
    real(wp) :: jmat(2,2) ! inverse of Jacobian

    real(wp) :: eps, jmat_pole(2,2), jmat_eps(2,2)

    eps = self % eps

    if ( eta(1) == 0.0_wp ) then

      jmat(:,:) = self % jmat_pole(:,:)

    else if ( 0.0_wp < eta(1) .and. eta(1) < eps ) then

      jmat_pole(:,:) = self % jmat_pole(:,:)

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
