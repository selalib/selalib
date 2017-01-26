module m_test_poisson_2d_polar_dirichlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_poisson_2d_polar_base, only: &
    c_test_poisson_2d_polar_base

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  implicit none

  public :: &
    t_test_poisson_2d_polar_dirichlet_quadratic, &
    t_test_poisson_2d_polar_dirichlet_cubic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  ! Dirichlet base class
  !-----------------------------------------------------------------------------
  type, extends(c_test_poisson_2d_polar_base), abstract :: c_test_dirichlet
  
    ! Boundary conditions (not overwritable)
    sll_int32, private :: bc_rmin = sll_p_dirichlet
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
    procedure :: get_rlim => dirichlet_get_rlim
    procedure :: get_bcs  => dirichlet_get_bcs

  end type c_test_dirichlet

  !-----------------------------------------------------------------------------
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = (r-rmax)(r-rmin)(a + b*cos(k(th-th0))
  !-----------------------------------------------------------------------------
  type, extends(c_test_dirichlet) :: t_test_poisson_2d_polar_dirichlet_quadratic

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => dirichlet_zero_error_phi_ex
    procedure :: phi_ex_diff1_r  => dirichlet_zero_error_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => dirichlet_zero_error_phi_ex_d2r
    procedure :: phi_ex_diff2_th => dirichlet_zero_error_phi_ex_d2th

  end type t_test_poisson_2d_polar_dirichlet_quadratic

  !-----------------------------------------------------------------------------
  ! Test-case with expected non-zero numerical error (cubic radial profile)
  ! \phi(r,th) = r(r-rmax)(r-rmin)(a + b*cos(k(th-th0))
  !-----------------------------------------------------------------------------
  type, extends(c_test_dirichlet) :: t_test_poisson_2d_polar_dirichlet_cubic

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => dirichlet_phi_ex
    procedure :: phi_ex_diff1_r  => dirichlet_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => dirichlet_phi_ex_d2r
    procedure :: phi_ex_diff2_th => dirichlet_phi_ex_d2th

  end type t_test_poisson_2d_polar_dirichlet_cubic

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

  !=============================================================================
  ! Test-case with expected zero numerical error (parabolic radial profile)
  ! \phi(r,th) = (r-rmax)(r-rmin)(a + b*cos(k(th-th0))
  !=============================================================================

  pure function dirichlet_zero_error_phi_ex( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_quadratic), intent(in) :: self
    sll_real64                                        , intent(in) :: r
    sll_real64                                        , intent(in) :: th
    sll_real64 :: val

    val = (r-self%rmax)*(r-self%rmin)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex

  !-----------------------------------------------------------------------------
  pure function dirichlet_zero_error_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_quadratic), intent(in) :: self
    sll_real64                                        , intent(in) :: r
    sll_real64                                        , intent(in) :: th
    sll_real64 :: val

    val = (2.0_f64*r-self%rmin-self%rmax)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex_d1r

  !-----------------------------------------------------------------------------
  pure function dirichlet_zero_error_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_quadratic), intent(in) :: self
    sll_real64                                        , intent(in) :: r
    sll_real64                                        , intent(in) :: th
    sll_real64 :: val

    val = 2.0_f64*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex_d2r

  !-----------------------------------------------------------------------------
  pure function dirichlet_zero_error_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_quadratic), intent(in) :: self
    sll_real64                                        , intent(in) :: r
    sll_real64                                        , intent(in) :: th
    sll_real64 :: val

    val = (r-self%rmax)*(r-self%rmin)*(-self%b*self%k**2*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex_d2th

  !=============================================================================
  ! Test-case with expected non-zero numerical error (cubic radial profile)
  ! \phi(r,th) = r(r-rmax)(r-rmin)(a + b*cos(k(th-th0))
  !=============================================================================

  pure function dirichlet_phi_ex( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_cubic), intent(in) :: self
    sll_real64                                    , intent(in) :: r
    sll_real64                                    , intent(in) :: th
    sll_real64 :: val

    val = r*(r-self%rmax)*(r-self%rmin)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex

  !-----------------------------------------------------------------------------
  pure function dirichlet_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_cubic), intent(in) :: self
    sll_real64                                    , intent(in) :: r
    sll_real64                                    , intent(in) :: th
    sll_real64 :: val

    val = (3.0_f64*r**2-2.0_f64*(self%rmin+self%rmax)*r+self%rmin*self%rmax) &
          *(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex_d1r

  !-----------------------------------------------------------------------------
  pure function dirichlet_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_cubic), intent(in) :: self
    sll_real64                                    , intent(in) :: r
    sll_real64                                    , intent(in) :: th
    sll_real64 :: val

    val = (6.0_f64*r-2.0_f64*(self%rmin+self%rmax)) &
          *(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex_d2r

  !-----------------------------------------------------------------------------
  pure function dirichlet_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_poisson_2d_polar_dirichlet_cubic), intent(in) :: self
    sll_real64                                    , intent(in) :: r
    sll_real64                                    , intent(in) :: th
    sll_real64 :: val

    val = r*(r-self%rmax)*(r-self%rmin)*(-self%b*self%k**2*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex_d2th

end module m_test_poisson_2d_polar_dirichlet
