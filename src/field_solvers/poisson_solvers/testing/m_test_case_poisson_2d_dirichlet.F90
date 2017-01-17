module m_test_case_poisson_2d_dirichlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_case_poisson_2d_base, only: &
    c_test_case_poisson_2d_polar

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  implicit none

  public :: &
    t_test_dirichlet_zero_error, &
    t_test_dirichlet, &
    sll_s_test_dirichlet_init

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(c_test_case_poisson_2d_polar) :: t_test_dirichlet_zero_error

    !---------------------------------------------------------------------------
    ! \phi(r,th) = (r-rmax)(r-rmin)(a + b*cos(k(th-th0))
    !---------------------------------------------------------------------------
    sll_real64 :: a   = 0.1_f64
    sll_real64 :: b   = 1.6_f64
    sll_real64 :: th0 = 0.3_f64
    sll_int32  :: k   = 3

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => dirichlet_zero_error_phi_ex
    procedure :: phi_ex_diff1_r  => dirichlet_zero_error_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => dirichlet_zero_error_phi_ex_d2r
    procedure :: phi_ex_diff2_th => dirichlet_zero_error_phi_ex_d2th

  end type t_test_dirichlet_zero_error

  type, extends(c_test_case_poisson_2d_polar) :: t_test_dirichlet

    !---------------------------------------------------------------------------
    ! \phi(r,th) = r(r-rmax)(r-rmin)(a + b*cos(k(th-th0))
    !---------------------------------------------------------------------------
    sll_real64 :: a   = 0.1_f64
    sll_real64 :: b   = 1.6_f64
    sll_real64 :: th0 = 0.3_f64
    sll_int32  :: k   = 3

  contains
    ! 2D manufactured solution
    procedure :: phi_ex          => dirichlet_phi_ex
    procedure :: phi_ex_diff1_r  => dirichlet_phi_ex_d1r
    procedure :: phi_ex_diff2_r  => dirichlet_phi_ex_d2r
    procedure :: phi_ex_diff2_th => dirichlet_phi_ex_d2th

  end type t_test_dirichlet

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine sll_s_test_dirichlet_init( self, rmin, rmax )
    class(c_test_case_poisson_2d_polar), intent(inout) :: self
    sll_real64, intent(in) :: rmin
    sll_real64, intent(in) :: rmax

    self%rmin = rmin
    self%rmax = rmax
    self%bc_rmin = sll_p_dirichlet
    self%bc_rmax = sll_p_dirichlet

  end subroutine sll_s_test_dirichlet_init

  !-----------------------------------------------------------------------------
  pure function dirichlet_zero_error_phi_ex( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (r-self%rmax)*(r-self%rmin)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex

  !-----------------------------------------------------------------------------
  pure function dirichlet_zero_error_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (2.0_f64*r-self%rmin-self%rmax)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex_d1r

  !-----------------------------------------------------------------------------
  pure function dirichlet_zero_error_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = 2.0_f64*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex_d2r

  !-----------------------------------------------------------------------------
  pure function dirichlet_zero_error_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (r-self%rmax)*(r-self%rmin)*(-self%b*self%k**2*cos( self%k*(th-self%th0) ))

  end function dirichlet_zero_error_phi_ex_d2th

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function dirichlet_phi_ex( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = r*(r-self%rmax)*(r-self%rmin)*(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex

  !-----------------------------------------------------------------------------
  pure function dirichlet_phi_ex_d1r( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (3.0_f64*r**2-2.0_f64*(self%rmin+self%rmax)*r+self%rmin*self%rmax) &
          *(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex_d1r

  !-----------------------------------------------------------------------------
  pure function dirichlet_phi_ex_d2r( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (6.0_f64*r-2.0_f64*(self%rmin+self%rmax)) &
          *(self%a+self%b*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex_d2r

  !-----------------------------------------------------------------------------
  pure function dirichlet_phi_ex_d2th( self, r, th ) result( val )
    class(t_test_dirichlet), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = r*(r-self%rmax)*(r-self%rmin)*(-self%b*self%k**2*cos( self%k*(th-self%th0) ))

  end function dirichlet_phi_ex_d2th

end module m_test_case_poisson_2d_dirichlet
