!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Test-cases with Neumann mode 0 BCs for quasi-neutrality equation in 2D polar coordinates.
!>
!> @details
!> This module contains:
!> - Neumann mode 0 base class with default parameters (all overwritable except for BCs)
!> - Test-case with parabolic radial profile:
!>   phi(r,theta) = a(r-rmax)(r-2rmin+rmax) + b(r-rmax)(r-rmin)cos(k(theta-theta_0))
!>
!> The numerical error expected on this test-case with a solver employing a
!> 2nd-order numerical method in the radial direction is zero (machine precision),
!> if k (Fourier mode) is <= ntheta/2.
!>

module m_test_qn_solver_2d_polar_neumann_mode0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_qn_solver_2d_polar_base, only: &
    c_test_qn_solver_2d_polar_base

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann_mode_0

  implicit none

  public :: &
    t_test_qn_solver_2d_polar_neumann_mode0_quadratic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  ! Neumann-mode-0 base class
  !-----------------------------------------------------------------------------
  type, extends(c_test_qn_solver_2d_polar_base), abstract :: c_test_neumann_mode0

    ! Boundary conditions (not overwritable)
    sll_int32, private :: bc_rmin = sll_p_neumann_mode_0
    sll_int32, private :: bc_rmax = sll_p_dirichlet

    ! Default parameters (may be overwritten by user)
    sll_real64 :: rmin = 1.0_f64
    sll_real64 :: rmax = 10.0_f64
    sll_real64 :: a    = 0.1_f64
    sll_real64 :: b    = 1.6_f64
    sll_real64 :: th0  = 0.3_f64
    sll_int32  :: k    = 3
    logical    :: adiabatic_electrons = .true.
    logical    :: use_zonal_flow = .true.
    sll_real64 :: epsilon_0 = 1.0_f64

  contains

    ! Get domain limits and boundary conditions
    procedure :: get_rlim       => neumann_mode0_get_rlim
    procedure :: get_bcs        => neumann_mode0_get_bcs
    procedure :: get_parameters => neumann_mode0_get_parameters

  end type c_test_neumann_mode0

  !-----------------------------------------------------------------------------
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = a(r-rmax)(r-2rmin+rmax) + b(r-rmax)(r-rmin)cos(k(th-th0))
  !-----------------------------------------------------------------------------
  type, extends(c_test_neumann_mode0) :: t_test_qn_solver_2d_polar_neumann_mode0_quadratic

  contains
    ! 1D input profiles
    procedure :: rho_m0          => f_test__rho_m0
    procedure :: rho_m0_diff1_r  => f_test__rho_m0_diff1_r
    procedure :: b_magn          => f_test__b_magn
    procedure :: b_magn_diff1_r  => f_test__b_magn_diff1_r
    procedure :: lambda          => f_test__lambda
    ! 2D manufactured solution
    procedure :: phi_ex          => f_test__phi_ex
    procedure :: phi_ex_avg_th   => f_test__phi_ex_avg_th
    procedure :: phi_ex_diff1_r  => f_test__phi_ex_diff1_r
    procedure :: phi_ex_diff2_r  => f_test__phi_ex_diff2_r
    procedure :: phi_ex_diff2_th => f_test__phi_ex_diff2_th

  end type t_test_qn_solver_2d_polar_neumann_mode0_quadratic

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

  !-----------------------------------------------------------------------------
  pure subroutine neumann_mode0_get_parameters( self, adiabatic_electrons, &
                                                use_zonal_flow, epsilon_0 )
    class(c_test_neumann_mode0), intent(in   ) :: self
    logical                    , intent(  out) :: adiabatic_electrons
    logical                    , intent(  out) :: use_zonal_flow
    sll_real64                 , intent(  out) :: epsilon_0

    adiabatic_electrons = self%adiabatic_electrons
    use_zonal_flow      = self%use_zonal_flow
    epsilon_0           = self%epsilon_0

  end subroutine neumann_mode0_get_parameters

  !=============================================================================
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = a(r-rmax)(r-2rmin+rmax) + b(r-rmax)(r-rmin)cos(k(th-th0))
  !=============================================================================

  !-----------------------------------------------------------------------------
  ! 1D PROFILES
  !-----------------------------------------------------------------------------

  pure function f_test__rho_m0( self, r ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64, intent(in) :: r
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax )
      val = 10.0_f64*(rmax-r)/(rmax-rmin) + 2.0_f64*(r-rmin)/(rmax-rmin)
    end associate

  end function f_test__rho_m0

  !-----------------------------------------------------------------------------
  pure function f_test__rho_m0_diff1_r( self, r ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64 :: val

    associate( rmin => self%rmin, rmax => self%rmax )
      val = -8.0_f64/(rmax-rmin)
    end associate

  end function f_test__rho_m0_diff1_r 

  !-----------------------------------------------------------------------------
  pure function f_test__b_magn( self, r ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64 :: val

    val = 1.0_f64

  end function f_test__b_magn

  !-----------------------------------------------------------------------------
  pure function f_test__b_magn_diff1_r( self, r ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64 :: val

    val = 0.0_f64

  end function f_test__b_magn_diff1_r

  !-----------------------------------------------------------------------------
  pure function f_test__lambda( self, r ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64 :: val

    associate( rmean => 0.5_f64*(self%rmax+self%rmin), &
               width => 0.1_f64*(self%rmax-self%rmin) )
      val = 0.01_f64 * (1.0_f64 + exp( -(r-rmean)**2/width**2 ))
    end associate

  end function f_test__lambda

  !-----------------------------------------------------------------------------
  ! 2D MANUFACTURED SOLUTION
  !-----------------------------------------------------------------------------

  pure function f_test__phi_ex( self, r, th ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64                                              , intent(in) :: th
    sll_real64 :: val

    val = self%a*(r-self%rmax)*(r-2.0_f64*self%rmin+self%rmax) &
          + self%b*(r-self%rmax)*(r-self%rmin)*cos( self%k*(th-self%th0) )

  end function f_test__phi_ex

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_avg_th( self, r, th ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64                                              , intent(in) :: th
    sll_real64 :: val

    if (self%k == 0) then
      val = self%a*(r-self%rmax)*(r-2.0_f64*self%rmin+self%rmax) &
            + self%b*(r-self%rmax)*(r-self%rmin)
    else
      val = self%a*(r-self%rmax)*(r-2.0_f64*self%rmin+self%rmax)
    end if

  end function f_test__phi_ex_avg_th

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_diff1_r( self, r, th ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64                                              , intent(in) :: th
    sll_real64 :: val

    val = self%a*2.0_f64*(r-self%rmin) &
          + self%b*(2.0_f64*r-self%rmin-self%rmax)*cos( self%k*(th-self%th0) )

  end function f_test__phi_ex_diff1_r

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_diff2_r( self, r, th ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64                                              , intent(in) :: th
    sll_real64 :: val

    val = self%a*2.0_f64 + self%b*2.0_f64*cos( self%k*(th-self%th0) )

  end function f_test__phi_ex_diff2_r

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_diff2_th( self, r, th ) result( val )
    class(t_test_qn_solver_2d_polar_neumann_mode0_quadratic), intent(in) :: self
    sll_real64                                              , intent(in) :: r
    sll_real64                                              , intent(in) :: th
    sll_real64 :: val

    val = -self%b*(r-self%rmax)*(r-self%rmin)*self%k**2*cos( self%k*(th-self%th0) )

  end function f_test__phi_ex_diff2_th

end module m_test_qn_solver_2d_polar_neumann_mode0
