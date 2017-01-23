module m_test_case_qn_2d_dirichlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_case_qn_2d_base, only: &
    c_test_case_qn_solver_2d_polar

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  implicit none

  public :: &
    t_test_dirichlet_zero_error

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  ! Dirichlet base class
  !-----------------------------------------------------------------------------
  type, extends(c_test_case_qn_solver_2d_polar), abstract :: c_test_dirichlet

    ! Boundary conditions (not overwritable)
    sll_int32, private :: bc_rmin = sll_p_dirichlet
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
    procedure :: get_rlim       => dirichlet_get_rlim
    procedure :: get_bcs        => dirichlet_get_bcs
    procedure :: get_parameters => dirichlet_get_parameters

  end type c_test_dirichlet

  !-----------------------------------------------------------------------------
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = (r-rmax)(r-rmin)(a + b*cos(k(th-th0))
  !-----------------------------------------------------------------------------
  type, extends(c_test_dirichlet) :: t_test_dirichlet_zero_error

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

  end type t_test_dirichlet_zero_error

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !=============================================================================
  ! Dirichlet base class
  !=============================================================================

  pure function dirichlet_get_rlim( self ) result( rlim )
    class(c_test_dirichlet), intent(in) :: self
    sll_real64 :: rlim(2)

    rlim(1) = self%rmin
    rlim(2) = self%rmax

  end function dirichlet_get_rlim

  !-----------------------------------------------------------------------------
  pure function dirichlet_get_bcs( self ) result( bcs )
    class(c_test_dirichlet), intent(in) :: self
    sll_int32 :: bcs(2)

    bcs(1) = self%bc_rmin
    bcs(2) = self%bc_rmax

  end function dirichlet_get_bcs

  !-----------------------------------------------------------------------------
  pure subroutine dirichlet_get_parameters( self, adiabatic_electrons, &
                                            use_zonal_flow, epsilon_0 )
    class(c_test_dirichlet), intent(in  ) :: self
    logical,    intent(out) :: adiabatic_electrons
    logical,    intent(out) :: use_zonal_flow
    sll_real64, intent(out) :: epsilon_0

    adiabatic_electrons = self%adiabatic_electrons
    use_zonal_flow      = self%use_zonal_flow
    epsilon_0           = self%epsilon_0

  end subroutine dirichlet_get_parameters

  !=============================================================================
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = (r-rmax)(r-rmin)(a + b*cos(k(th-th0))
  !=============================================================================

  !-----------------------------------------------------------------------------
  ! 1D PROFILES
  !-----------------------------------------------------------------------------

  pure function f_test__rho_m0( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    associate( rmin => self%rmin, rmax => self%rmax )
      val = 10.0_f64*(rmax-r)/(rmax-rmin) + 2.0_f64*(r-rmin)/(rmax-rmin)
    end associate
  end function f_test__rho_m0

  !-----------------------------------------------------------------------------
  pure function f_test__rho_m0_diff1_r( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    associate( rmin => self%rmin, rmax => self%rmax )
      val = -8.0_f64/(rmax-rmin)
    end associate
  end function f_test__rho_m0_diff1_r 

  !-----------------------------------------------------------------------------
  pure function f_test__b_magn( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    val = 1.0_f64
  end function f_test__b_magn

  !-----------------------------------------------------------------------------
  pure function f_test__b_magn_diff1_r( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    val = 0.0_f64
  end function f_test__b_magn_diff1_r

  !-----------------------------------------------------------------------------
  pure function f_test__lambda( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
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
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (r-self%rmax)*(r-self%rmin)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function f_test__phi_ex

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_avg_th( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    if (self%k == 0) then
      val = (r-self%rmax)*(r-self%rmin)*(self%a+self%b)
    else
      val = (r-self%rmax)*(r-self%rmin)*self%a
    end if

  end function f_test__phi_ex_avg_th

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_diff1_r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (2.0_f64*r-self%rmin-self%rmax)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function f_test__phi_ex_diff1_r

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_diff2_r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = 2.0_f64*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function f_test__phi_ex_diff2_r

  !-----------------------------------------------------------------------------
  pure function f_test__phi_ex_diff2_th( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (r-self%rmax)*(r-self%rmin)*(-self%b*self%k**2*cos( self%k*(th-self%th0) ))

  end function f_test__phi_ex_diff2_th

end module m_test_case_qn_2d_dirichlet
