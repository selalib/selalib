module m_test_case_qn_2d_dirichlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use m_test_case_qn_2d_base, only: c_test_case_qn_solver_2d_polar

  implicit none

  public :: t_test_dirichlet_zero_error

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(c_test_case_qn_solver_2d_polar) :: t_test_dirichlet_zero_error

    ! Angular profile is a+cos(k*(th-th0))
    sll_real64, private :: a   = 0.1_f64
    sll_int32 , private :: k   = 3
    sll_real64, private :: th0 = 0.3_f64

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

  ! 1D profiles
  pure function f_test__rho_m0( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    associate( rmin => self%rmin, rmax => self%rmax )
      val = 10.0_f64*(rmax-r)/(rmax-rmin) + 2.0_f64*(r-rmin)/(rmax-rmin)
    end associate
  end function f_test__rho_m0


  pure function f_test__rho_m0_diff1_r( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    associate( rmin => self%rmin, rmax => self%rmax )
      val = -8.0_f64/(rmax-rmin)
    end associate
  end function f_test__rho_m0_diff1_r 


  pure function f_test__b_magn( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    val = 1.0_f64
  end function f_test__b_magn


  pure function f_test__b_magn_diff1_r( self, r ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64 :: val
    val = 0.0_f64
  end function f_test__b_magn_diff1_r


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
  ! 2D manufactured solution
  pure function f_test__phi_ex( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (self%rmax-r) * (r-self%rmin) * (self%a + cos( self%k*(th-self%th0) ))

  end function f_test__phi_ex


  pure function f_test__phi_ex_avg_th( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (self%rmax-r) * (r-self%rmin) * self%a

  end function f_test__phi_ex_avg_th


  pure function f_test__phi_ex_diff1_r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = -2.0_f64*r + self%rmin + self%rmax
    val = val * (self%a + cos( self%k*(th-self%th0) ))

  end function f_test__phi_ex_diff1_r


  pure function f_test__phi_ex_diff2_r( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = -2.0_f64 * (self%a + cos( self%k*(th-self%th0) ))

  end function f_test__phi_ex_diff2_r


  pure function f_test__phi_ex_diff2_th( self, r, th ) result( val )
    class(t_test_dirichlet_zero_error), intent(in) :: self
    sll_real64                        , intent(in) :: r
    sll_real64                        , intent(in) :: th
    sll_real64 :: val

    val = (self%rmax-r) * (r-self%rmin)
    val = val * real( -self%k**2, f64 ) * cos( self%k*(th-self%th0) )

  end function f_test__phi_ex_diff2_th


end module m_test_case_qn_2d_dirichlet
