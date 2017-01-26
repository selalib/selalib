!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Test-cases with Neumann mode 0 BCs for Poisson's equation in 2D polar coordinates.
!>
!> @details
!> This module contains:
!> - Neumann mode 0 base class with default parameters (all overwritable except for BCs)
!> - Test-case with parabolic radial profile:
!>   phi(r,theta) = a(r-rmax)(r-2rmin+rmax) + b(r-rmax)(r-rmin)cos(k(theta-theta_0))
!>
!> The numerical error expected on this test-case with a solver employing a
!> 2nd-order numerical method in the radial direction is zero (machine precision).
!>

module m_test_poisson_2d_polar_neumann_mode0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_poisson_2d_polar_base, only: &
    c_test_poisson_2d_polar_base

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann_mode_0

  implicit none

  public :: &
    t_test_poisson_2d_polar_neumann_mode0_quadratic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  ! Neumann-mode-0 base class
  !-----------------------------------------------------------------------------
  type, extends(c_test_poisson_2d_polar_base), abstract :: c_test_neumann_mode0

    ! Boundary conditions (non overwritable)
    sll_int32, private :: bc_rmin = sll_p_neumann_mode_0
    sll_int32, private :: bc_rmax = sll_p_dirichlet

    ! Default parameters (may be overwritten by user)
    sll_real64 :: rmin = 1.0_f64
    sll_real64 :: rmax = 2.0_f64
    sll_real64 :: a    = 0.1_f64
    sll_real64 :: b    = 1.6_f64
    sll_real64 :: th0  = 0.3_f64
    sll_int32  :: k    = 3

  contains

    ! Get domain limits and boundary conditions
    procedure :: get_rlim => neumann_mode0_get_rlim
    procedure :: get_bcs  => neumann_mode0_get_bcs

  end type c_test_neumann_mode0

  !-----------------------------------------------------------------------------
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = a(r-rmax)(r-2rmin+rmax) + b(r-rmax)(r-rmin)cos(k(th-th0))
  !-----------------------------------------------------------------------------
  type, extends(c_test_neumann_mode0) :: t_test_poisson_2d_polar_neumann_mode0_quadratic

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => neumann_mode0_zero_error_phi_ex
    procedure :: phi_ex_diff1_r  => neumann_mode0_zero_error_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => neumann_mode0_zero_error_phi_ex_d2r
    procedure :: phi_ex_diff2_th => neumann_mode0_zero_error_phi_ex_d2th

  end type t_test_poisson_2d_polar_neumann_mode0_quadratic

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !=============================================================================
  ! Neumann-mode-0 base class
  !=============================================================================

  pure function neumann_mode0_get_rlim( self ) result( rlim )
    class(c_test_neumann_mode0), intent(in) :: self
    sll_real64 :: rlim(2)

    rlim(1) = self%rmin
    rlim(2) = self%rmax

  end function neumann_mode0_get_rlim

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_get_bcs( self ) result( bcs )
    class(c_test_neumann_mode0), intent(in) :: self
    sll_int32 :: bcs(2)

    bcs(1) = self%bc_rmin
    bcs(2) = self%bc_rmax

  end function neumann_mode0_get_bcs

  !=============================================================================
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = a(r-rmax)(r-2rmin+rmax) + b(r-rmax)(r-rmin)cos(k(th-th0))
  !=============================================================================

  pure function neumann_mode0_zero_error_phi_ex( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                            , intent(in) :: r
    sll_real64                                            , intent(in) :: th
    sll_real64 :: val

    val = self%a*(r-self%rmax)*(r-2.0_f64*self%rmin+self%rmax) &
          + self%b*(r-self%rmax)*(r-self%rmin)*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_zero_error_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                            , intent(in) :: r
    sll_real64                                            , intent(in) :: th
    sll_real64 :: val

    val = self%a*2.0_f64*(r-self%rmin) &
          + self%b*(2.0_f64*r-self%rmin-self%rmax)*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex_d1r

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_zero_error_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                            , intent(in) :: r
    sll_real64                                            , intent(in) :: th
    sll_real64 :: val

    val = self%a*2.0_f64 + self%b*2.0_f64*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex_d2r

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_zero_error_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                            , intent(in) :: r
    sll_real64                                            , intent(in) :: th
    sll_real64 :: val

    val = -self%b*(r-self%rmax)*(r-self%rmin)*self%k**2*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex_d2th

end module m_test_poisson_2d_polar_neumann_mode0
