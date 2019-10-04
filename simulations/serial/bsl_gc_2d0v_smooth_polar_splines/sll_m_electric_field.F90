module sll_m_electric_field
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_pi

  use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

  use sll_m_spline_2d, only: sll_t_spline_2d

  implicit none

  public :: sll_t_electric_field

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_electric_field

    ! Can be overwritten at initialization
    real(wp) :: eps = 1.0e-12_wp

    type(sll_t_singular_mapping_discrete), pointer :: mapping_discrete => null()
    type(sll_t_spline_2d)                , pointer :: spline_2d_phi    => null()

  contains

    procedure :: init => s_electric_field__init
    procedure :: eval => s_electric_field__eval
    procedure :: free => s_electric_field__free

  end type sll_t_electric_field

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_electric_field__init( self, mapping_discrete, spline_2d_phi, eps )
    class(sll_t_electric_field)                   , intent(inout) :: self
    type (sll_t_singular_mapping_discrete), target, intent(in   ) :: mapping_discrete
    type (sll_t_spline_2d)                , target, intent(in   ) :: spline_2d_phi
    real(wp), optional                            , intent(in   ) :: eps

    self % mapping_discrete => mapping_discrete

    self % spline_2d_phi => spline_2d_phi

    if ( present( eps ) ) self % eps = eps

  end subroutine s_electric_field__init

  !-----------------------------------------------------------------------------
  ! Implements strategy described in section 5.4 of https://doi.org/10.1016/j.jcp.2019.108889
  !-----------------------------------------------------------------------------
  SLL_PURE function s_electric_field__eval( self, eta ) result( electric_field )
    class(sll_t_electric_field), intent(in) :: self
    real(wp)                   , intent(in) :: eta(2)
    real(wp) :: electric_field(2)

    real(wp) :: jdet, d1, d2, d3, d4, d5, d6, th1, th2, eps, ef_0(2), ef_eps(2), jmat(2,2)

    eps = self % eps

    if ( eta(1) == 0.0_wp ) then

      th1 = 0.0_wp
      th2 = 0.5_wp * sll_p_pi

      d1 = self % spline_2d_phi % eval_deriv_x1( eta(1), th1 ) ! dphi/ds(0,theta_1)
      d2 = self % spline_2d_phi % eval_deriv_x1( eta(1), th2 ) ! dphi/ds(0,theta_2)

      jmat = self % mapping_discrete % jmat( (/ eta(1), th1 /) )
 
      d3 = jmat(1,1) ! dx/ds(0,theta_1)
      d4 = jmat(2,1) ! dy/ds(0,theta_1)

      jmat = self % mapping_discrete % jmat( (/ eta(1), th2 /) )
 
      d5 = jmat(1,1) ! dx/ds(0,theta_2)
      d6 = jmat(2,1) ! dy/ds(0,theta_2)

      electric_field(1) = - ( d4*d2 - d1*d6 ) / ( d4*d5 - d3*d6 )
      electric_field(2) = - ( d1*d5 - d3*d2 ) / ( d4*d5 - d3*d6 )

    ! Implements last equation of section 5.4 of https://doi.org/10.1016/j.jcp.2019.108889
    else if ( 0.0_wp < eta(1) .and. eta(1) < eps ) then

      th1 = 0.0_wp
      th2 = 0.5_wp * sll_p_pi

      d1 = self % spline_2d_phi % eval_deriv_x1( 0.0_wp, th1 ) ! dphi/ds(0,theta_1)
      d2 = self % spline_2d_phi % eval_deriv_x1( 0.0_wp, th2 ) ! dphi/ds(0,theta_2)

      jmat = self % mapping_discrete % jmat( (/ 0.0_wp, th1 /) )
 
      d3 = jmat(1,1) ! dx/ds(0,theta_1)
      d4 = jmat(2,1) ! dy/ds(0,theta_1)

      jmat = self % mapping_discrete % jmat( (/ 0.0_wp, th2 /) )
 
      d5 = jmat(1,1) ! dx/ds(0,theta_2)
      d6 = jmat(2,1) ! dy/ds(0,theta_2)

      ! E(0,theta)
      ef_0(1) = - ( d4*d2 - d1*d6 ) / ( d4*d5 - d3*d6 )
      ef_0(2) = - ( d1*d5 - d3*d2 ) / ( d4*d5 - d3*d6 )

      jmat = self % mapping_discrete % jmat( (/ eps, eta(2) /) )
      jdet = self % mapping_discrete % jdet( (/ eps, eta(2) /) )

      ! E(eps,theta): J^(-T) times 'logical' gradient
      ef_eps(1) = - (   jmat(2,2) * self % spline_2d_phi % eval_deriv_x1( eps, eta(2) ) &
                      - jmat(2,1) * self % spline_2d_phi % eval_deriv_x2( eps, eta(2) ) ) / jdet
      ef_eps(2) = - ( - jmat(1,2) * self % spline_2d_phi % eval_deriv_x1( eps, eta(2) ) &
                      + jmat(1,1) * self % spline_2d_phi % eval_deriv_x2( eps, eta(2) ) ) / jdet

      ! Linear interpolation between 0 and eps
      electric_field(1) = (1.0_wp-eta(1)/eps)*ef_0(1) + eta(1)/eps*ef_eps(1)
      electric_field(2) = (1.0_wp-eta(1)/eps)*ef_0(2) + eta(1)/eps*ef_eps(2)

    ! Implements equation (29) of https://doi.org/10.1016/j.jcp.2019.108889
    else if ( eps <= eta(1) ) then

      jmat = self % mapping_discrete % jmat( eta )
      jdet = self % mapping_discrete % jdet( eta )

      ! J^(-T) times 'logical' gradient
      electric_field(1) = - (   jmat(2,2) * self % spline_2d_phi % eval_deriv_x1( eta(1), eta(2) ) &
                              - jmat(2,1) * self % spline_2d_phi % eval_deriv_x2( eta(1), eta(2) ) ) / jdet
      electric_field(2) = - ( - jmat(1,2) * self % spline_2d_phi % eval_deriv_x1( eta(1), eta(2) ) &
                              + jmat(1,1) * self % spline_2d_phi % eval_deriv_x2( eta(1), eta(2) ) ) / jdet

    end if

  end function s_electric_field__eval

  !-----------------------------------------------------------------------------
  subroutine s_electric_field__free( self )
    class(sll_t_electric_field), intent(inout) :: self

    nullify( self % mapping_discrete )

    nullify( self % spline_2d_phi )

  end subroutine s_electric_field__free

end module sll_m_electric_field
