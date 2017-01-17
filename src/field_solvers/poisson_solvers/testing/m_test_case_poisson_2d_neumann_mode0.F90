module m_test_case_poisson_2d_neumann_mode0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_case_poisson_2d_base, only: &
    c_test_case_poisson_2d_polar

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann_mode_0

  implicit none

  public :: &
    t_test_neumann_mode0_zero_error, &
    sll_s_test_neumann_mode0_init

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(c_test_case_poisson_2d_polar) :: t_test_neumann_mode0_zero_error

    !---------------------------------------------------------------------------
    ! \phi(r,th) = a(r-rmax)(r-2rmin+rmax) + b(r-rmax)(r-rmin)cos(k(th-th0))
    !---------------------------------------------------------------------------
    sll_real64 :: a   = 0.1_f64
    sll_real64 :: b   = 1.6_f64
    sll_real64 :: th0 = 0.3_f64
    sll_int32  :: k   = 3

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => neumann_mode0_zero_error_phi_ex
    procedure :: phi_ex_diff1_r  => neumann_mode0_zero_error_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => neumann_mode0_zero_error_phi_ex_d2r
    procedure :: phi_ex_diff2_th => neumann_mode0_zero_error_phi_ex_d2th

  end type t_test_neumann_mode0_zero_error

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine sll_s_test_neumann_mode0_init( self, rmin, rmax )
    class(c_test_case_poisson_2d_polar), intent(inout) :: self
    sll_real64, intent(in) :: rmin
    sll_real64, intent(in) :: rmax

    self%rmin = rmin
    self%rmax = rmax
    self%bc_rmin = sll_p_neumann_mode_0
    self%bc_rmax = sll_p_dirichlet

  end subroutine sll_s_test_neumann_mode0_init

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_zero_error_phi_ex( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = self%a*(r-self%rmax)*(r-2.0_f64*self%rmin+self%rmax) &
          + self%b*(r-self%rmax)*(r-self%rmin)*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_zero_error_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = self%a*2.0_f64*(r-self%rmin) &
          + self%b*(2.0_f64*r-self%rmin-self%rmax)*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex_d1r

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_zero_error_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = self%a*2.0_f64 + self%b*2.0_f64*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex_d2r

  !-----------------------------------------------------------------------------
  pure function neumann_mode0_zero_error_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_neumann_mode0_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = -self%b*(r-self%rmax)*(r-self%rmin)*self%k**2*cos( self%k*(th-self%th0) )

  end function neumann_mode0_zero_error_phi_ex_d2th

end module m_test_case_poisson_2d_neumann_mode0
